##########################################################
#	File: Config.py
#	Author: Jose Perez <josegperez@mail.com>
#	Version: Model v5
#
#	Modify this file to your heart's content
#	Here we have constants used by my code
#	I put them all in one place so they can be easy to modify
#
##########################################################

# How much CD28 is needed to become active
CD28_THRESHOLD = 2

# How much time a cell will wait before resetting (mcs)
WAIT_TIME = 10

# How old the cells are when born (volume)
INITIAL_AGE = 15.0

# How old the cells will get at every step
STEP_AGE = 1.0
# How old the cells need to be to choose between apoptosis, division, and quiescence
DECISION_AGE = 50.0

# Probabilities of what the cell will do stochastically once it reaches a certain age/size
# Must add up to 1.0
PROB_APOPTOSIS = 0.1
PROB_DIVISION = 0.4
PROB_QUIESCENCE = 0.5

# Probabilities of which ligand will be lost when losing a ligand
PROB_LOST_CD80 = 0.5
PROB_LOST_CD86 = 0.5

# Probabilities of receptors binding to CD80
# Must add up to 1.0
PROB_CTLA4_BIND_CD80 = .9524
PROB_CD28_BIND_CD80 = .0476
WEIGHTS_CD80 = [PROB_CTLA4_BIND_CD80, PROB_CD28_BIND_CD80]

# Probabilities of receptors binding to CD86
# Must add up to 1.0
PROB_CTLA4_BIND_CD86 = .8837
PROB_CD28_BIND_CD86 = .1163
WEIGHTS_CD86 = [PROB_CTLA4_BIND_CD86, PROB_CD28_BIND_CD86]