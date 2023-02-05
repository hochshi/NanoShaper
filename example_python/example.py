import NanoShaper

PyINFO = ' <<PyINFO>>'
PyERROR = ' <<PyERROR>>'

# Example of usage of NanoShaper from Python

##################### init and conf file ################
conf = 'conf32.prm'
print "\n",PyINFO,"NanoShaper version "+ NanoShaper.VERSION +" Python binding starting"

# set it to false if atom info is absent (see file 1rszH.pdb.xyzr)
# in this mode on log you get the info
hasAtomInfo = False
# set it to true if atom info is present (see file 1pyt_plus.xyzr)
# in this mode you get info in ProShape format
#hasAtomInfo = True

# consistency check of conf file and conf loading
cf = NanoShaper.init(conf)

scale = cf.readFloat("Grid_scale")
perfil = cf.readFloat("Grid_perfil")
molFile = cf.readString("XYZR_FileName")
buildStatus = False

if (cf.readString("Build_status_map").lower()=='true'):
	buildStatus=True

########################################################

# normalMode and pocketMode are the two library functions of NanoShaper
if (cf.readString("Operative_Mode") == 'normal'):
	grid = NanoShaper.DelPhiShared(scale,perfil,molFile,False,buildStatus,False);
	surf = NanoShaper.surfaceFactory().create(cf,grid);
	NanoShaper.normalMode(surf,grid)
	del grid
	del surf
elif (cf.readString("Operative_Mode") == 'pockets'):
	NanoShaper.pocketMode(hasAtomInfo,cf);
else:
	print '\n',PyERROR,"Unknown operative mode"

NanoShaper.dispose(cf)
print ""