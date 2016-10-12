from pymol import cmd

def setAquariaColors():
	cmd.set_color("aquaria_sheet", [0.96, 0.78, 0.09])
	cmd.set_color("aquaria_helix", [0.39, 0.568, 0.714])
	cmd.set_color("aquaria_coil", [0.435, 0.619, 0.321])

def colorByAquaria(selection='all'):
	setAquariaColors()
	cmd.color("aquaria_helix", "ss H in " + selection)
	cmd.color("aquaria_sheet", "ss S in " + selection)
	cmd.color("aquaria_coil", "ss L+'' in " + selection)

cmd.extend('colorByAquaria', colorByAquaria)
cmd.extend('setAquariaColors', setAquariaColors)