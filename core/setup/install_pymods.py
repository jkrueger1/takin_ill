#
# installs Takin's Magdyn and instrument data loader python modules
# @author Tobias Weber <tweber@ill.fr>
# @date 16-may-2024
# @license GPLv2
#

import os
import site
import shutil


magdyn     = [ "pymods/magdyn.py", "pymods/_magdyn_py.so" ]
instr      = [ "pymods/instr.py", "pymods/_instr_py.so" ]


sitepath = site.getusersitepackages()
print("Using user's site package directory: \"%s\"." % sitepath)

if not os.path.isdir(sitepath):
	print("Site packages directory does not yet exist, creating it.")
	os.makedirs(sitepath)


print("Installing Takin/Magdyn...")
for file in magdyn:
	print("\t\"%s\" -> \"%s\"" % (file, sitepath))
	shutil.copy(file, sitepath)

print("Installing Takin/Instr...")
for file in instr:
	print("\t\"%s\" -> \"%s\"" % (file, sitepath))
	shutil.copy(file, sitepath)
