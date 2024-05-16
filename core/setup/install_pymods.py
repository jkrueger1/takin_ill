#
# installs Takin's Magdyn and instrument data loader python modules
# @author Tobias Weber <tweber@ill.fr>
# @date 16-may-2024
# @license GPLv2
#

import os
import site
import shutil


create_links = True

moddir     = "pymods/"
magdyn     = [ "magdyn.py", "_magdyn_py.so" ]
instr      = [ "instr.py", "_instr_py.so" ]


sitepath = site.getusersitepackages()
print("Using user's site package directory: \"%s\"." % sitepath)

if not os.path.isdir(sitepath):
	print("Site packages directory does not yet exist, creating it.")
	os.makedirs(sitepath)


def install_mod(files, msg):
	print(msg)

	for file in files:
		file_src = moddir + file
		file_dst = sitepath + "/" + file

		print("\t\"%s\" -> \"%s\"" % (file_src, file_dst))

		os.remove(file_dst)

		if create_links:
			os.symlink(os.path.abspath(file_src), file_dst)
		else:
			shutil.copy(file_src, sitepath)


install_mod(magdyn, "Installing Takin/Magdyn...")
install_mod(instr, "Installing Takin/Instr...")
