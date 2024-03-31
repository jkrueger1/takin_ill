#
# plots data files
# @author Tobias Weber <tweber@ill.fr>
# @date 1-apr-2024
# @license see 'LICENSE' file
#
# ----------------------------------------------------------------------------
# mag-core (part of the Takin software suite)
# Copyright (C) 2018-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
#                          Grenoble, France).
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

import sys
import math
import instr


#
# normalise detector (y) by monitor (m) counts
# y_new = y / m
# dy_new = 1/m dy - y/m^2 dm
#
def norm_counts_to_mon(y, dy, m, dm):
	val = y / m;
	err = math.sqrt(math.pow(dy/m, 2.) + math.pow(dm*y/(m*m), 2.))

	return [val, err]


# get energy transfer from ki and kf
def get_E(ki, kf):
	#E_to_k2 = 2.*co.neutron_mass/hbar_in_meVs**2. / co.elementary_charge*1000. * 1e-20
	E_to_k2 = 0.482596406464  # calculated with scipy, using the formula above

	return (ki**2. - kf**2.) / E_to_k2


def load_data(datfile, mergefiles = []):
	print("Loading \"%s\"." % (datfile))
	dat = instr.FileInstrBaseD.LoadInstr(datfile)
	if dat == None:
		return

	for mergefile in mergefiles:
		tomerge = instr.FileInstrBaseD.LoadInstr(mergefile)
		if tomerge != None:
			print("Merging with \"%s\"." % (mergefile))
			dat.MergeWith(tomerge)


	hs = []
	ks = []
	ls = []
	Es = []
	Is = []
	Is_err = []

	cntcol = dat.GetCol(dat.GetCountVar())
	moncol = dat.GetCol(dat.GetMonVar())

	# iterate scan points
	for point_idx in range(cntcol.size()):
		(h, k, l, ki, kf) = dat.GetScanHKLKiKf(point_idx)
		E = get_E(ki, kf)

		hs.append(h)
		ks.append(k)
		ls.append(l)
		Es.append(E)

		counts = cntcol[point_idx]
		counts_err = math.sqrt(counts)
		mon = moncol[point_idx]
		mon_err = math.sqrt(mon)

		if counts == 0:
			counts_err = 1.
		if mon == 0:
			mon_err = 1.

		[I, I_err] = norm_counts_to_mon(counts, counts_err, mon, mon_err)
		Is.append(I)
		Is_err.append(I_err)

	return (hs, ks, ls, Es, Is, Is_err)


# load scan files
if len(sys.argv) < 2:
	print("No scan file given.")
	exit(-1)

scanfile = sys.argv[1]
mergefiles = sys.argv[2:]
[hs, ks, ls, Es, Is, Is_err] = load_data(scanfile, mergefiles)


# plot the scan data
import matplotlib.pyplot as plt

plt.figure()
plt.xlabel("E (meV)")
plt.ylabel("S (a.u.)")
plt.errorbar(Es, Is, Is_err, marker='o', capsize=4, ls="none")
plt.tight_layout()
plt.show()
