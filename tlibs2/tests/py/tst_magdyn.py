#
# tlibs2 python interface test
# @author Tobias Weber <tweber@ill.fr>
# @date 12-oct-2023
# @license GPLv3, see 'LICENSE' file
#
# ----------------------------------------------------------------------------
# tlibs
# Copyright (C) 2017-2023  Tobias WEBER (Institut Laue-Langevin (ILL),
#                          Grenoble, France).
# Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
#                          (TUM), Garching, Germany).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# ----------------------------------------------------------------------------
#

import numpy
import tl2_magdyn


# files
modelfile = "../../../data/demos/magnon_Cu2OSeO3/model.magdyn"
dispfile = "disp.dat"


# create the magdyn object
mag = tl2_magdyn.MagDynD()


# load the model file
print("Loading {}...".format(modelfile))
if mag.Load(modelfile):
	print("Loaded {}.".format(modelfile))
else:
	print("Failed loading {}.".format(modelfile))
	exit(-1)


# directly calculate a dispersion and write it to a file
#print("Goldstone Energy: {:.4f} meV".format(mag.GetGoldstoneEnergy()))
print("\nSaving dispersion to {}...".format(dispfile))
mag.SaveDispersion(dispfile,  0, 0, 0.5,  1, 1, 0.5,  128)


# manually calculate the same dispersion and output it to the console
print("\nManual calculation of a dispersion...")
print("{:>15} {:>15} {:>15} {:>15} {:>15}".format("h", "k", "l", "E", "S(Q,E)"))
for h in numpy.linspace(0, 1, 128):
	k = h
	l = 0.5
	for S in mag.GetEnergies(h, k, l, False):
		print("{:15.4f} {:15.4f} {:15.4f} {:15.4f} {:15.4g}".format(h, k, l, S.E, S.weight))
