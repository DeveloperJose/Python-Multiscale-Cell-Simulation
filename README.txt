APC & TCell Model v5
By Jose Perez <josegperez@mail.com>

Description:
	This is a simple model created in CompuCell3D that simulates the interaction between Antigen Presenting Cells and T Cells.
	For more information about the project please read my abstract and poster.
	This project was presented at the 2016 Annual Biomedical Engineering Society Conference.

Requirements:
	CompuCell3D w/ twedit++
	Basic knowledge of XML
	Some knowledge of Python
	Some knowledge of the CC3D documentation

How To Use:
	Open TCELL_APC_Model.cc3d using twedit++

Notes:
	I have tried to document every single line in all of my files.
	You might need to go over some Python tutorials to understand the comments.
	If you would like a more general overview of the project please read my abstract or my poster.
	You can email me if you have any questions.


	This model doesn't have the SBML file for the CTLA-4 recycling.
	By looking at one of the examples (deltaNotch with SBML) you can easily create a simple one.
	BioNetGen also generates SBML files but they are SBML v3 and CC3D only supports up to SBML v2 files.