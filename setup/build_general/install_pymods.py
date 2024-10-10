#
# installs takin's python modules
# @author Tobias Weber <tweber@ill.fr>
# @date 16-may-2024
# @license GPLv2
#

import os
import sys
import site
import shutil


magdyn     = [ "magdyn.py", "_magdyn_py.so" ]
instr      = [ "instr.py",  "_instr_py.so" ]
bz         = [ "bzcalc.py", "_bz_py.so" ]


# get args
create_links = True

for arg in sys.argv[1:]:
	if arg == "-l":
		create_links = True
	elif arg == "-c":
		create_links = False


sitepath = site.getusersitepackages()
print("Using user's site package directory: \"%s\"." % sitepath)
if create_links:
	print("Installing by creating symlinks.")
else:
	print("Installing by copying files.")


modpath = sys.path[0] + "/pymods/"
if not os.path.isdir(modpath):
	modpath = sys.path[0] + "/../pymods/"

if not os.path.isdir(sitepath):
	print("Site packages directory does not yet exist, creating it.")
	os.makedirs(sitepath)


def install_mod(files, msg):
	print(msg)

	for file in files:
		file_src = modpath + file
		file_dst = sitepath + "/" + file

		print("\t\"%s\" -> \"%s\"." % (file_src, file_dst))

		if os.path.isfile(file_dst):
			os.remove(file_dst)

		if create_links:
			os.symlink(os.path.abspath(file_src), file_dst)
		else:
			shutil.copy(file_src, sitepath)


install_mod(magdyn, "\nInstalling Takin/Magdyn...")
install_mod(instr, "\nInstalling Takin/Instr...")
install_mod(bz, "\nInstalling Takin/BZ...")
print()
