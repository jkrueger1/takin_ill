#
# magdyn py interface demo -- plotting multiple dispersion branches
# @author Tobias Weber <tweber@ill.fr>
# @date 18-sep-2024
# @license GPLv3, see 'LICENSE' file
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

import numpy
import magdyn


# options
print_dispersion       = False  # write dispersion to console
only_positive_energies = True   # ignore magnon annihilation?

# dispersion branches to plot
dispersion = [
	numpy.array([ 0., 0., 0.5 ]),
	numpy.array([ 1., 1., 0.5 ]),
	numpy.array([ 1., 0.5, 1. ]),
]

# Q components to use for plotting branches
dispersion_plot_indices = [ 0, 2 ]

# number of Qs to calculate on a dispersion branch
num_Q_points = 256

# weight scaling and clamp factors
S_scale = 64.
S_min   = 1.
S_max   = 500.

# files
modelfile = "model.magdyn"


# create the magdyn object
mag = magdyn.MagDyn()


# load the model file
print("Loading {}...".format(modelfile))
if not mag.Load(modelfile):
	print("Failed loading {}.".format(modelfile))
	exit(-1)


# minimum energy
print("Energy minimum at Q = (000): {:.4f} meV".format(mag.CalcMinimumEnergy()))
print("Ground state energy: {:.4f} meV".format(mag.CalcGroundStateEnergy()))


# number of dispersion braches
num_branches = len(dispersion) - 1

# calculate the dispersion branches
print("\nCalculating %d dispersion branches..." % num_branches)
if print_dispersion:
	print("{:>15} {:>15} {:>15} {:>15} {:>15}".format("h", "k", "l", "E", "S(Q,E)"))

# data for all branches
data = []

for disp_idx in range(num_branches):
	hkl_start = dispersion[disp_idx]
	hkl_end = dispersion[disp_idx + 1]

	print("[%d/%d] Calculating %s  ->  %s branch..." %
	   (disp_idx + 1, num_branches, hkl_start, hkl_end))

	data_h = []
	data_k = []
	data_l = []
	data_E = []
	data_S = []

	# calculate a branch
	for hkl in numpy.linspace(hkl_start, hkl_end, num_Q_points):
		for S in mag.CalcEnergies(hkl[0], hkl[1], hkl[2], False):
			if only_positive_energies and S.E < 0.:
				continue

			weight = S.weight * S_scale
			if weight < S_min:
				weight = S_min
			elif weight > S_max:
				weight = S_max

			data_h.append(hkl[0])
			data_k.append(hkl[1])
			data_l.append(hkl[2])
			data_E.append(S.E)
			data_S.append(weight)

			if print_dispersion:
				print("{:15.4f} {:15.4f} {:15.4f} {:15.4f} {:15.4g}".format(
					hkl[0], hkl[1], hkl[2], S.E, S.weight))

	data.append([ disp_idx, data_h, data_k, data_l, data_E, data_S ])


# plot the results
print("\nPlotting dispersion branches...")

import matplotlib.pyplot as plot

(plt, axes) = plot.subplots(nrows = 1, ncols = num_branches, sharey = True)

for [ disp_idx, data_h, data_k, data_l, data_E, data_S ] in data:
	# dispersion branch start and end points
	b1 = dispersion[disp_idx]
	b2 = dispersion[disp_idx + 1]

	plot_idx = dispersion_plot_indices[disp_idx]
	if plot_idx == 0:
		data_x = data_h
		label = "(h %.4g %.4g) -> (h %.4g %.4g)" % (b1[1], b1[2], b2[1], b2[2])
		label_x = "h (rlu)"
	elif plot_idx == 1:
		data_x = data_k
		label = "(%.4g k %.4g) -> (%.4g k %.4g)" % (b1[0], b1[2], b2[0], b2[2])
		label_x = "k (rlu)"
	elif plot_idx == 2:
		data_x = data_l
		label = "(%.4g %.4g l) -> (%.4g %.4g l)" % (b1[0], b1[1], b2[0], b2[1])
		label_x = "l (rlu)"

	# labels
	axes[disp_idx].set_xlabel(label + ", " + label_x)
	if disp_idx == 0:
		axes[disp_idx].set_ylabel("E (meV)")

	# plot the dispersion branch
	axes[disp_idx].scatter(data_x, data_E, marker = '.', s = data_S)

plt.tight_layout()
plt.subplots_adjust(wspace = 0)
#plot.savefig("model.pdf")
plot.show()
