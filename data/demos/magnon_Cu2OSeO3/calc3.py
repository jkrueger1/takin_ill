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
import concurrent.futures as fut
import matplotlib.pyplot as plot


# -----------------------------------------------------------------------------
# options
print_dispersion       = False  # write dispersion to console
only_positive_energies = True   # ignore magnon annihilation?

# number of worker threads
max_threads = 8

# dispersion branches to plot
dispersion = [
	numpy.array([ 0., 0., 0.5 ]),
	numpy.array([ 1., 1., 0.5 ]),
	numpy.array([ 1., 0.5, 1. ]),
]

# number of Qs to calculate on a dispersion branch
num_Q_points = 256

# weight scaling and clamp factors
S_scale      = 64.
S_clamp_min  = 1.
S_clamp_max  = 500.
S_filter_min = -1.    # don't filter

# files
modelfile = "model.magdyn"
plotfile  = ""
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# create the magdyn object
mag = magdyn.MagDyn()


#
# calculate the energies and weights for a Q point
#
def calc_Es(h, k, l):
	Es_dicts = mag.CalcEnergies(h, k, l, False)

	Es = []
	for Es_dict in Es_dicts:
		Es.append(( h, k, l, Es_dict.E, Es_dict.weight ))

	return Es


# load the model file
print("Loading {}...".format(modelfile))
if not mag.Load(modelfile):
	print("Failed loading {}.".format(modelfile))
	exit(-1)


# minimum energy
print("Energy minimum at Q = (000): {:.4f} meV".format(mag.CalcMinimumEnergy()))
print("Ground state energy: {:.4f} meV".format(mag.CalcGroundStateEnergy()))
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# number of dispersion braches
num_branches = len(dispersion) - 1

# Q components to use for plotting branches
dispersion_plot_indices = []

# calculate the dispersion branches
print("\nCalculating %d dispersion branches..." % num_branches)
if print_dispersion:
	print("{:>15} {:>15} {:>15} {:>15} {:>15}".format("h", "k", "l", "E", "S(Q,E)"))

# data for all branches
data = []

for branch_idx in range(num_branches):
	hkl_start = dispersion[branch_idx]
	hkl_end = dispersion[branch_idx + 1]

	print("[%d/%d] Calculating %s  ->  %s branch..." %
	   (branch_idx + 1, num_branches, hkl_start, hkl_end))

	# find scan axis
	Q_diff = [
		numpy.abs(hkl_start[0] - hkl_end[0]),
		numpy.abs(hkl_start[1] - hkl_end[1]),
		numpy.abs(hkl_start[2] - hkl_end[2]) ]

	axis_idx = 0
	if Q_diff[1] > Q_diff[axis_idx]:
		axis_idx = 1
	elif Q_diff[2] > Q_diff[axis_idx]:
		axis_idx = 2
	dispersion_plot_indices.append(axis_idx)

	data_h = []
	data_k = []
	data_l = []
	data_E = []
	data_S = []

	# calculate a branch
	Es_futures = []
	with fut.ProcessPoolExecutor(max_workers = max_threads) as exe:
		# submit tasks
		for hkl in numpy.linspace(hkl_start, hkl_end, num_Q_points):
			Es_futures.append(exe.submit(calc_Es, hkl[0], hkl[1], hkl[2]))

		# get results from tasks
		for Es_future in Es_futures:
			for (h, k, l, E, weight) in Es_future.result():
				if only_positive_energies and E < 0.:
					continue

				weight = weight * S_scale

				if weight > S_filter_min:
					if weight < S_clamp_min:
						weight = S_clamp_min
					elif weight > S_clamp_max:
						weight = S_clamp_max

					data_h.append(h)
					data_k.append(k)
					data_l.append(l)
					data_E.append(E)
					data_S.append(weight)

				if print_dispersion:
					print("{:15.4f} {:15.4f} {:15.4f} {:15.4f} {:15.4g}".format(
						h, k, l, E, weight))

	data.append([ branch_idx, data_h, data_k, data_l, data_E, data_S ])
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# plot the results
print("\nPlotting dispersion branches...")

(plt, axes) = plot.subplots(nrows = 1, ncols = num_branches, sharey = True)

for [ branch_idx, data_h, data_k, data_l, data_E, data_S ] in data:
	# dispersion branch start and end points
	b1 = dispersion[branch_idx]
	b2 = dispersion[branch_idx + 1]

	plot_idx = dispersion_plot_indices[branch_idx]
	if plot_idx == 0:
		data_x = data_h
		detailed_label = "(h %.4g %.4g) -> (h %.4g %.4g)" % (b1[1], b1[2], b2[1], b2[2])
		label_x = "h (rlu)"
	elif plot_idx == 1:
		data_x = data_k
		detailed_label = "(%.4g k %.4g) -> (%.4g k %.4g)" % (b1[0], b1[2], b2[0], b2[2])
		label_x = "k (rlu)"
	elif plot_idx == 2:
		data_x = data_l
		detailed_label = "(%.4g %.4g l) -> (%.4g %.4g l)" % (b1[0], b1[1], b2[0], b2[1])
		label_x = "l (rlu)"

	# ticks and labels
	axes[branch_idx].set_xlim(data_x[0], data_x[-1])
	if branch_idx == 0:
		axes[branch_idx].set_ylabel("E (meV)")
		tick_labels = [
			"(%.4g %.4g %.4g)" % (b1[0], b1[1], b1[2]),
			"(%.4g %.4g %.4g)" % (b2[0], b2[1], b2[2])]
	else:
		tick_labels = ["", "(%.4g %.4g %.4g)" % (b2[0], b2[1], b2[2])]
	axes[branch_idx].set_xticks([data_x[0], data_x[-1]], labels = tick_labels)

	if branch_idx == num_branches / 2 - 1:
		axes[branch_idx].set_xlabel("Q (rlu)")

	# detailed label
	#axes[branch_idx].set_xlabel(detailed_label + ", " + label_x)

	# plot the dispersion branch
	axes[branch_idx].scatter(data_x, data_E, marker = '.', s = data_S)

	# flip x axis
	#if b2[plot_idx] < b1[plot_idx]:
	#	axes[branch_idx].invert_xaxis()

plt.tight_layout()
plt.subplots_adjust(wspace = 0)

if plotfile != "":
	plot.savefig(plotfile)
plot.show()
# -----------------------------------------------------------------------------