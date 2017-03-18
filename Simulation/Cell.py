##########################################################
#	File: Cell.py
#	Author: Jose Perez <josegperez@mail.com>
#	Version: Model v5
#
#	This is where the magic happens.
#	This is the code for the interactions between APCs and TCells
#	In this model version all interactions (intracellular and extracellular) happen inside CC3D
#	The SBML model for intracellular activity is commented out as I am fixing it
#	
##########################################################
import Data
import Config

# This is just to have a global constant for our dictionaries
class CC3DKey(object):
    DATA_KEY = 'DATA'

# This is to be able to change the type of cell that CC3D displays easily from inside this file.
# This is set-up in the MainSteppable so the values aren't actually -1
class CC3DType(object):
    MEDIUM = 0
    
    APC = -1
    
    TREG_INACTIVE = -1
    TREG_ACTIVE = -1
    TREG_ANERGIC = -1
    
    TCONV_INACTIVE = -1
    TCONV_ACTIVE = -1
    TCONV_ANERGIC = -1

# The state of a TCell
class State(object):
    ACTIVE = 0
    INACTIVE = 1
    ANERGIC = 2
    
	# Used for APC-TCell ligand and receptor binding.
	# When Peptide-MHC and TCR bind they call others to the binding area
	# This is part of a process called TCell co-activation.
	# Please see Kaur14PhD (http://etheses.bham.ac.uk/4903/8/Kaur14PhD.pdf) for more information about co-activation.
    AWAITING_COACTIVATION = 3

# This is our parent class for our cells
# At first it was going to be used to make sure no duplicate code was used in my code
# However at the moment it's really just a nice visualization of the cell methods.
class CellData(object):    
    def __init__(self, cc3d_cell):
		# We store the CC3D cell type inside so we can change it.
        self.cc3d_cell = cc3d_cell
        
    def interact_with_apc(self, apc, mcs):
        pass
        
    def interact_with_tcell(self, tcell, mcs):
        pass
        
    def reset(self, initialize=False):
        pass

		
# Antigen Presenting Cell
class APC(CellData):       
    def __init__(self, cc3d_cell):
		# Call our parent's constructor
        CellData.__init__(self, cc3d_cell)

        # Set initial/default quantities for APCs
        self.initial_PEPTIDEMHC = 10
        self.initial_CD80 = 15
        self.initial_CD86 = 15
        
        # Set our APC to default
        self.reset(initialize=True)
        
		# Increase our count of APCs for the plots
        Data.TOTAL_APC += 1
 
    def reset(self, initialize=False):
		# This method resets APCs
		# It's meant to be used every time an APC stops interacting with a TCell
		# It unbinds everything basically
		# Initialize is there to let me know if we are creating the APC or if we just finished an interaction
		# True = We are creating the cell
		# False = We are resetting the cell
		
        if not initialize:
            # ----=== Global Data ===---- #
			# To update the plots
            Data.TOTAL_AMOUNT_PEPTIDEMHC -= self.total_PEPTIDEMHC
            Data.TOTAL_AMOUNT_CD80 -= self.total_CD80
            Data.TOTAL_AMOUNT_CD86 -= self.total_CD86
            
        self.total_PEPTIDEMHC = self.initial_PEPTIDEMHC
        self.total_CD80 = self.initial_CD80
        self.total_CD86 = self.initial_CD86
        
        # ----=== Global Data ===---- #
        Data.TOTAL_AMOUNT_PEPTIDEMHC += self.total_PEPTIDEMHC
        Data.TOTAL_AMOUNT_CD80 += self.total_CD80
        Data.TOTAL_AMOUNT_CD86 += self.total_CD86
    
    def remove_ligand(self):
		# Remove a random ligand from this APC
		# We don't actually remove ligands, we just simulate doing it
		# You can check the plot of lost ligands to see how many ligands get would be lost hypothetically
        import random
        
		# Let's see which ligands are unbound and store them here
        ligands_available = []
        ligand = ''
        
		# If we have CD80, then we can remove some
        if self.initial_CD80 > 0:
            ligands_available.append('CD80')
            
		# If we have CD86, then we can remove some
        if self.initial_CD86 > 0:
            ligands_available.append('CD86')
        
        # We can't remove ligands if there are none
        if not ligands_available:
            return
        
        # If we only have one ligand then we'll remove that one
        if len(ligands_available) == 1:
            ligand = ligands_available[0]
        
        # Choose randomly which ligand will be taken away if there is more than one
		# We just throw a dice. We could use weights if we wanted to (see MitosisSteppable)
        if not ligand:  
            dice = random.randrange(0, 2)
            ligand = 'CD80' if dice == 0 else 'CD86'
            
        self._internal_remove(ligand)
        
    def _internal_remove(self, ligand):
        if ligand == 'CD80':
			# To actually remove CD80 uncomment the line below
			# self.initial_CD80 -= 1
            # ----=== Global Data ===---- #  
            Data.TOTAL_LOST_CD80 += 1
        
        elif ligand == 'CD86':
			# To actually remove CD86 uncomment the line below
			# self.initial_CD86 -= 1
            # ----=== Global Data ===---- #  
            Data.TOTAL_LOST_CD86 += 1

# TCell class
# Most of the magic happens here
class TCell(CellData):  
	# Types of TCells
    TREG = 0
    TCONV = 1

    # Ligands (CD80, CD86)
    # Receptors (CTLA-4, CD28)
        
    def __init__(self, cc3d_cell, type, state = State.INACTIVE):
        CellData.__init__(self, cc3d_cell)
        
		# Set our type to the specified type
        self.type = type
		# Set our state to the specified State
		# Default is State.INACTIVE
        self.state = state
        
		# This is used to simulate if the TCell is internalizing CTLA-4
		# If True we are internalizing CTLA-4, otherwise False
        self.int = False
        
		# Set-up our initial concentrations
        self.reset(initialize=True)
        # ----=== Global Data for plots ===---- #
        Data.TOTAL_TCELLS += 1
        if self.type == self.TREG:
            Data.TOTAL_TREG_INACTIVE += 1 if state == State.INACTIVE else 0
            Data.TOTAL_TREG_ACTIVE += 1 if state == State.ACTIVE else 0
            Data.TOTAL_TREG_ANERGIC += 1 if state == State.ANERGIC else 0
        else:
            Data.TOTAL_TCONV_INACTIVE += 1 if state == State.INACTIVE else 0
            Data.TOTAL_TCONV_ACTIVE += 1 if state == State.ACTIVE else 0
            Data.TOTAL_TCONV_ANERGIC += 1 if state == State.ANERGIC else 0
            
    def reset(self, initialize=False):
		# What's the ID of the APC we are bound to?
		# -1 if we aren't bound to an APC
        self.bound_to_id = -1
		# If we are bound to an APC this will contain all the info of that APC
        self.bound_to = 0
        
		# What was the last MCS we were bound?
		self.bound_last_mcs = -1
		# How long have we been bound to an APC?
        self.bound_time = 0
		# How long have we been unbound to an APC?
        self.unbound_time = 0
       
		# How much CD28 is bound to an APC?
        self.bound_CD28 = 0
        
        if not initialize:
            # ----=== Global Data ===---- # 
            Data.TOTAL_AMOUNT_TCR -= self.total_TCR
            Data.TOTAL_AMOUNT_CD28 -= self.total_CD28
            Data.TOTAL_AMOUNT_EXTERNAL_CTLA4 -= self.total_external_CTLA4
            Data.TOTAL_AMOUNT_INTERNAL_CTLA4 -= self.total_internal_CTLA4
        
        # Initial variables for tregs
        if self.type == self.TREG:
            self.total_TCR = 50
            self.total_CD28 = 25            
            self.total_external_CTLA4 = 10
            self.total_internal_CTLA4 = 0
        # Initial variables for everything else (tconv)
        else:
            self.total_TCR = 50
            self.total_CD28 = 25            
            self.total_external_CTLA4 = 0
            self.total_internal_CTLA4 = 0 # Added during activation    

        # ----=== Global Data ===---- # 
        Data.TOTAL_AMOUNT_TCR += self.total_TCR
        Data.TOTAL_AMOUNT_CD28 += self.total_CD28
        Data.TOTAL_AMOUNT_EXTERNAL_CTLA4 += self.total_external_CTLA4
        Data.TOTAL_AMOUNT_INTERNAL_CTLA4 += self.total_internal_CTLA4
        
        
    def log(self, message):
		# I use this method to debug errors in the code
		# It prints the cell type plus the id plus a specified message
        typeFormat = 'TREG' if self.type == self.TREG else 'TCONV'
        print '[' + typeFormat + ' Cell #' + str(self.cc3d_cell.id) + '] ' + str(message)
        
    def contact_with_friend(self, apc, mcs):
		# This method is called when we touch an APC
        # Check if we aren't bound to make a friend
		# AKA binding to that APC
        if self.bound_to_id == -1:
			# Bind to that APC!
            self.bound_to_id = apc.cc3d_cell.id
            self.bound_to = apc
            return True
        # Check we are talking to our same friend
		# We are going to ignore other APCs and just interact with our "friends"
        elif self.bound_to_id != apc.cc3d_cell.id:
            # Start missing interaction timer
			# This will start counting how long it's been since we last interacted with our "friend"
            if self.bound_last_mcs == -1:
                self.bound_last_mcs = mcs
            else:
				# Our timer is running. Calculate how much time has passed.
				# deltaTime = currentTime - lastRecordedTime
                time = mcs - self.bound_last_mcs
				
				# Time has passed, so store it
                if time >= 1:
                    self.bound_last_mcs = mcs                   
                    self.unbound_time += time                   
                    
            # Reset if we haven't interacted with our friend
			# We do this if more time has passed than we are patient enough to wait
            if self.unbound_time >= Config.WAIT_TIME:
                #self.log('Was bound for about ' + str(self.bound_time))
				# Unbind from APC
                self.bound_time = 0
                self.bound_to.reset()
                self.reset()
                
                # Ligand removal
				# Perhaps the TCell took a ligand from the APC
                if self.type == self.TREG:
                    apc.remove_ligand()  
                
            return False
        
        else:
            #self.log('Bound for about ' + str(self.bound_time))
            return True  
            
    def interact_with_apc(self, apc, mcs):  
		# If the TCell is inactive
		# Bind the TCR
        if self.state == State.INACTIVE:
            self.bind_tcr(apc, mcs)
        
		# TCR was bound
		# Co-activation is required
        elif self.state == State.AWAITING_COACTIVATION:
            self.select_interaction(apc, mcs)

		# If the TCell is active, randomly internalize some CTLA-4
		# This requires more testing
		# Also you can replace this with the CTLA-4 recyling SBML
        elif self.state == State.ACTIVE:
			# If we are internalizing CTLA-4
            if self.int:
                self.total_internal_CTLA4 += 1
                Data.TOTAL_AMOUNT_INTERNAL_CTLA4 += 1
                
                self.total_external_CTLA4 -= 1
                Data.TOTAL_AMOUNT_EXTERNAL_CTLA4 -= 1
            # Otherwise we are externalizing CTLA-4
			else:
                self.total_internal_CTLA4 -= 1
                Data.TOTAL_AMOUNT_INTERNAL_CTLA4 -= 1
                
                self.total_external_CTLA4 += 1
                Data.TOTAL_AMOUNT_EXTERNAL_CTLA4 += 1
            
			# Toggle between internalizing and externalizing randomly
            if mcs % 10 == 0:
                self.int = not self.int
    
    def bind_tcr(self, apc, mcs):
        if not self.contact_with_friend(apc, mcs):
            return
            
        # We are binding to an APC
		# Reset our unbound time
        self.unbound_time = 0
        
		# Set-up our bound time starting from the current MCS
        time = mcs - self.bound_last_mcs
        if time >= 1:
            self.bound_last_mcs = mcs
            self.bound_time += time
        
        # Bind TCR and Peptide-MHC
		# If we have both available only
        if self.total_TCR > 0 and apc.total_PEPTIDEMHC > 0:
				# We require co-activation now
                self.state = State.AWAITING_COACTIVATION
				# Bind TCR and MHC
                self.total_TCR -= 1
                apc.total_PEPTIDEMHC -= 1
				# Try to bind some ligands and receptors
                self.select_interaction(apc, mcs)  
                
                # ----=== Global Data ===---- #
				# Update our plots
                Data.TOTAL_AMOUNT_TCR -= 1
                Data.TOTAL_AMOUNT_PEPTIDEMHC -= 1
    
    def select_interaction(self, apc, mcs):
        if not self.contact_with_friend(apc, mcs):
            return
        
		# In order to bind ligands and receptors we
		# have to first figure out which ones we still
		# have available
        ligands_available = []
        receptors_available = []
        
        ligand = ''
        receptor = ''
        
		# If we have CD86 then it can bind
        if apc.total_CD86 > 0:
            ligands_available.append('CD86')
        
		# If we have CD80 then it can bind
        if apc.total_CD80 > 0:
            ligands_available.append('CD80')
        
        if not ligands_available:
            # The cell needed co-stimulation but didn't receive it
            # Change T-Cell into Anergic T-Cell
			# This is based on Kaur14PhD's description of co-activation
            #self.log("Required co-stimulation but didn't receive it. Transforming into anergic.")
            if self.state == State.AWAITING_COACTIVATION:
				self.cc3d_cell.type = CC3DType.TREG_ANERGIC if self.type == self.TREG else CC3DType.TCONV_ANERGIC
				# Set the state to ANERGIC
                self.state = State.ANERGIC
                
				# Kill it
                self.cc3d_cell.targetVolume = 0
                self.cc3d_cell.lambdaVolume = 0
                
				# Update our plots
                Data.TOTAL_TCELLS -= 1
                Data.TOTAL_TREG_INACTIVE -= 1 if self.type == self.TREG else 0
                Data.TOTAL_TCONV_INACTIVE -= 1 if self.type == self.TCONV else 0
                Data.TOTAL_TREG_ANERGIC += 1 if self.type == self.TREG else 0
                Data.TOTAL_TCONV_ANERGIC += 1 if self.type == self.TCONV else 0
            return   
    
        # Select a random ligand
        from numpy.random import choice
        ligand = choice(ligands_available)    
        
        # NOTE: CLTA-4 MUST be appended first for the weights
        if self.total_external_CTLA4 > 0:
            receptors_available.append('CTLA-4')      
  
        if self.total_CD28 > 0:
            receptors_available.append('CD28')
        
        # No receptors available
        if not receptors_available:
            return
        
        # Only 1 receptor available, pick it
        if len(receptors_available) == 1:
            receptor = receptors_available[0]
        
        # If there was only 1 receptor...
        if receptor:
			# Then we select it
            self.match_with_apc(apc, ligand, receptor, mcs)
        else:
			# Otherwise randomly select based on weights
			# The weights are based on affinities from the Kaur14PhD paper
            weights = Config.WEIGHTS_CD80 if ligand == 'CD80' else Config.WEIGHTS_CD86
            receptor = choice(receptors_available, p=weights)
			
            self.match_with_apc(apc, ligand, receptor, mcs)        
    
    def match_with_apc(self, apc, ligand, receptor, mcs):
        # If we have CD28....
		if receptor == 'CD28' and self.total_CD28 > 0:
			# This is for our plots
            Data.TOTAL_ENGAGED_CD28 += 1
			# Count how much is bound for our activation threshold
            self.bound_CD28 += 1
            
			# Decrease our CD28 current
            self.total_CD28 -= 1
			# Update the plots
            Data.TOTAL_AMOUNT_CD28 -= 1
            
			# Decrease the quantities from the APC depending on what we are binding to
			# Either CD80 or CD86
            apc.total_CD80 -= 1 if ligand == 'CD80' else 0
            apc.total_CD86 -= 1 if ligand == 'CD86' else 0

			# Update our plots
            Data.TOTAL_AMOUNT_CD80 -= 1 if ligand == 'CD80' else 0
            Data.TOTAL_AMOUNT_CD86 -= 1 if ligand == 'CD86' else 0
            
            # Activate T-Cell if they pass the CD28 threshold
            if self.bound_CD28 > Config.CD28_THRESHOLD:
				# Change the CC3D simulation type to active
                self.cc3d_cell.type = CC3DType.TREG_ACTIVE if self.type == self.TREG else CC3DType.TCONV_ACTIVE
				# Change our internal state to active
                self.state = State.ACTIVE
                
				# Update our plots
                Data.TOTAL_TREG_INACTIVE -= 1 if self.type == self.TREG else 0
                Data.TOTAL_TCONV_INACTIVE -= 1 if self.type == self.TCONV else 0
                
                Data.TOTAL_TREG_ACTIVE += 1 if self.type == self.TREG else 0
                Data.TOTAL_TCONV_ACTIVE += 1 if self.type == self.TCONV else 0
                
				# Add CTLA-4 to TCONV once it becomes active
				if self.type == self.TCONV:
					self.total_external_CTLA4 += 1
					Data.TOTAL_AMOUNT_EXTERNAL_CTLA4 += 1
					self.total_internal_CTLA4 += 1
					Data.TOTAL_AMOUNT_INTERNAL_CTLA4 += 1
               
		# CTLA-4 binding
        elif receptor == 'CTLA-4' and self.total_external_CTLA4 > 0:
			# Update our plot
            Data.TOTAL_ENGAGED_EXTERNAL_CTLA4 += 1
            
			# Decrease external ctla-4 count
            self.total_external_CTLA4 -= 1
			# Update CTLA-4 plot
            Data.TOTAL_AMOUNT_EXTERNAL_CTLA4 -= 1
            
			# Decrease CD80 or CD86 depending on which one we are binding to
			# We decrease by 2 (based on figures of Kaur14PhD paper)
            apc.total_CD80 -= 2 if ligand == 'CD80' else 0
            apc.total_CD86 -= 2 if ligand == 'CD86' else 0

			# Update plots
            Data.TOTAL_AMOUNT_CD80 -= 2 if ligand == 'CD80' else 0
            Data.TOTAL_AMOUNT_CD86 -= 2 if ligand == 'CD86' else 0
            
			# Update our plots
			# This is so we don't go into the negatives
            if apc.total_CD80 < 0:
                Data.TOTAL_AMOUNT_CD80 += (-1 * apc.total_CD80)              
                apc.total_CD80 = 0
                
            if apc.total_CD86 < 0:
                Data.TOTAL_AMOUNT_CD80 += (-1 * apc.total_CD86)                 
                apc.total_CD86 = 0