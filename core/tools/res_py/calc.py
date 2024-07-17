#!/usr/bin/env python3
#
# calculates the instrumental resolution
#
# @author Tobias Weber <tweber@ill.fr>
# @date feb-2015, oct-2019, jul-2024
# @license GPLv2
#
# ----------------------------------------------------------------------------
# Takin (inelastic neutron scattering software package)
# Copyright (C) 2017-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
#                          Grenoble, France).
# Copyright (C) 2013-2017  Tobias WEBER (Technische Universitaet Muenchen
#                          (TUM), Garching, Germany).
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; version 2 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
# ----------------------------------------------------------------------------
#

# requires numpy version >= 1.10
import numpy as np
import helpers
import params_in20

np.set_printoptions(floatmode = "fixed",  precision = 4)


# -----------------------------------------------------------------------------
#
# default resolution parameters
# NOTE: not all parameters are used by all calculation backends
#
params = {
    # options
    "verbose" : True,

    # resolution method, "eck", "pop", or "cn"
    "reso_method" : "pop",

    # scattering triangle
    "ki" : 1.4,
    "kf" : 1.4,
    "E"  : helpers.get_E(1.4, 1.4),
    "Q"  : 1.777,

    # d spacings
    "mono_xtal_d" : 3.355,
    "ana_xtal_d"  : 3.355,

     # scattering senses
    "mono_sense"   : -1.,
    "sample_sense" :  1.,
    "ana_sense"    : -1.,
    "mirror_Qperp" : False,

    # distances
    "dist_src_mono"    :  10. * helpers.cm2A,
    "dist_mono_sample" : 200. * helpers.cm2A,
    "dist_sample_ana"  : 115. * helpers.cm2A,
    "dist_ana_det"     :  85. * helpers.cm2A,

    # shapes
    "src_shape"    : "rectangular",  # "rectangular" or "circular"
    "sample_shape" : "cylindrical",  # "cuboid" or "cylindrical"
    "det_shape"    : "rectangular",  # "rectangular" or "circular"

    # component sizes
    "src_w"    : 6.   * helpers.cm2A,
    "src_h"    : 12.  * helpers.cm2A,
    "mono_d"   : 0.15 * helpers.cm2A,
    "mono_w"   : 12.  * helpers.cm2A,
    "mono_h"   : 8.   * helpers.cm2A,
    "sample_d" : 1.   * helpers.cm2A,
    "sample_w" : 1.   * helpers.cm2A,
    "sample_h" : 1.   * helpers.cm2A,
    "ana_d"    : 0.3  * helpers.cm2A,
    "ana_w"    : 12.  * helpers.cm2A,
    "ana_h"    : 8.   * helpers.cm2A,
    "det_w"    : 1.5  * helpers.cm2A,
    "det_h"    : 5.   * helpers.cm2A,

    # horizontal collimation
    "coll_h_pre_mono"    : 30. * helpers.min2rad,
    "coll_h_pre_sample"  : 30. * helpers.min2rad,
    "coll_h_post_sample" : 30. * helpers.min2rad,
    "coll_h_post_ana"    : 30. * helpers.min2rad,

    # vertical collimation
    "coll_v_pre_mono"    : 30. * helpers.min2rad,
    "coll_v_pre_sample"  : 30. * helpers.min2rad,
    "coll_v_post_sample" : 30. * helpers.min2rad,
    "coll_v_post_ana"    : 30. * helpers.min2rad,

    # horizontal focusing
    "mono_curv_h" : 0.,
    "ana_curv_h"  : 0.,
    "mono_is_curved_h" : False,
    "ana_is_curved_h"  : False,
    "mono_is_optimally_curved_h" : False,
    "ana_is_optimally_curved_h"  : False,
    "mono_curv_h_formula" : None,
    "ana_curv_h_formula" : None,

    # vertical focusing
    "mono_curv_v" : 0.,
    "ana_curv_v"  : 0.,
    "mono_is_curved_v" : False,
    "ana_is_curved_v"  : False,
    "mono_is_optimally_curved_v" : False,
    "ana_is_optimally_curved_v"  : False,
    "mono_curv_v_formula" : None,
    "ana_curv_v_formula" : None,

    # guide before monochromator
    "use_guide"   : True,
    "guide_div_h" : 15. * helpers.min2rad,
    "guide_div_v" : 15. * helpers.min2rad,

    # horizontal mosaics
    "mono_mosaic"   : 45. * helpers.min2rad,
    "sample_mosaic" : 30. * helpers.min2rad,
    "ana_mosaic"    : 45. * helpers.min2rad,

    # vertical mosaics
    "mono_mosaic_v"   : 45. * helpers.min2rad,
    "sample_mosaic_v" : 30. * helpers.min2rad,
    "ana_mosaic_v"    : 45. * helpers.min2rad,

    # calculate R0 factor (not needed if only the ellipses are to be plotted)
    "calc_R0" : True,

    # crystal reflectivities; TODO, so far always 1
    "dmono_refl" : 1.,
    "dana_effic" : 1.,

    # off-center scattering
    # WARNING: while this is calculated, it is not yet considered in the ellipse plots
    "pos_x" : 0. * helpers.cm2A,
    "pos_y" : 0. * helpers.cm2A,
    "pos_z" : 0. * helpers.cm2A,

    # vertical scattering in kf, keep "False" for normal TAS
    "kf_vert" : False,
}
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# get command-line arguments
# -----------------------------------------------------------------------------
import argparse
argparser = argparse.ArgumentParser(
    description = "Calculates the resolution ellipsoid of a TAS instrument.")

argparser.add_argument("-i", "--instr", default = None, type = str,
    help = "set the parameters to a pre-defined instrument (in20fc)")
argparser.add_argument("--silent", action = "store_true",
    help = "disable output")
argparser.add_argument("-p", "--plot", action = "store_true",
    help = "plot results")
argparser.add_argument("-m", "--reso_method", default = None, type = str,
    help = "resolution method to use (cn/pop/eck)")
argparser.add_argument("--kf_vert", action = "store_true",
    help = "scatter vertically in the kf axis (only for Eckold-Sobolev method)")
argparser.add_argument("--ki", default = None, type = float,
    help = "incoming wavenumber")
argparser.add_argument("--kf", default = None, type = float,
    help = "outgoing wavenumber")
argparser.add_argument("-E", default = None, type = float,
    help = "energy transfer")
argparser.add_argument("-Q", default = None, type = float,
    help = "momentum transfer")
argparser.add_argument("--mono_sense", default = None, type = float,
    help = "monochromator scattering sense")
argparser.add_argument("--sample_sense", default = None, type = float,
    help = "sample scattering sense")
argparser.add_argument("--ana_sense", default = None, type = float,
    help = "analyser scattering sense")
argparser.add_argument("--mono_curv_v", default = None, type = float,
    help = "monochromator vertical curvature")
argparser.add_argument("--mono_curv_h", default = None, type = float,
    help = "monochromator horizontal curvature")
argparser.add_argument("--ana_curv_v", default = None, type = float,
    help = "analyser vertical curvature")
argparser.add_argument("--ana_curv_h", default = None, type = float,
    help = "analyser horizontal curvature")
argparser.add_argument("--mono_curv_v_opt", action = "store_true",
    help = "set optimal monochromator vertical curvature")
argparser.add_argument("--mono_curv_h_opt", action = "store_true",
    help = "set optimal monochromator horizontal curvature")
argparser.add_argument("--ana_curv_v_opt", action = "store_true",
    help = "set optimal analyser vertical curvature")
argparser.add_argument("--ana_curv_h_opt", action = "store_true",
    help = "set optimal analyser horizontal curvature")

parsedargs = argparser.parse_args()

if parsedargs.instr != None:
    if parsedargs.instr == "in20fc":
        print("Loaded IN20/Flatcone default parameters.\n")
        params = params_in20.params_fc

# get parsed command-line arguments
show_plots = parsedargs.plot
params["verbose"] = not parsedargs.silent

if parsedargs.reso_method != None:
    params["reso_method"] = parsedargs.reso_method

params["kf_vert"] = parsedargs.kf_vert

if parsedargs.ki != None:
    params["ki"] = parsedargs.ki
if parsedargs.kf != None:
    params["kf"] = parsedargs.kf

if parsedargs.Q != None:
    params["Q"] = parsedargs.Q

if parsedargs.mono_sense != None:
    params["mono_sense"] = parsedargs.mono_sense
if parsedargs.sample_sense != None:
    params["sample_sense"] = parsedargs.sample_sense
if parsedargs.ana_sense != None:
    params["ana_sense"] = parsedargs.ana_sense

# ensure that the senses are either +1 or -1
if params["mono_sense"] > 0.:
    params["mono_sense"] = 1.
else:
    params["mono_sense"] = -1.
if params["sample_sense"] > 0.:
    params["sample_sense"] = 1.
else:
    params["sample_sense"] = -1.
if params["ana_sense"] > 0.:
    params["ana_sense"] = 1.
else:
    params["ana_sense"] = -1.

# set instrument position
if parsedargs.E != None:
    # calculate ki from E and kf
    params["E"] = parsedargs.E
    params["ki"] = helpers.get_ki(params["E"], params["kf"])
else:
    # calculate E from ki and kf
    params["E"] = helpers.get_E(params["ki"], params["kf"])

# set fixed curvatures
if parsedargs.mono_curv_v != None:
    params["mono_curv_v"] = parsedargs.mono_curv_v
    params["mono_is_curved_v"] = True
if parsedargs.mono_curv_h != None:
    params["mono_curv_h"] = parsedargs.mono_curv_h
    params["mono_is_curved_h"] = True
if parsedargs.ana_curv_v != None:
    params["ana_curv_v"] = parsedargs.ana_curv_v
    params["ana_is_curved_v"] = True
if parsedargs.ana_curv_h != None:
    params["ana_curv_h"] = parsedargs.ana_curv_h
    params["ana_is_curved_h"] = True

# set optimal curvatures
if parsedargs.mono_curv_v_opt:
    params["mono_is_optimally_curved_v"] = True
    params["mono_is_curved_v"] = True
if parsedargs.mono_curv_h_opt:
    params["mono_is_optimally_curved_h"] = True
    params["mono_is_curved_h"] = True
if parsedargs.ana_curv_v_opt:
    params["ana_is_optimally_curved_v"] = True
    params["ana_is_curved_v"] = True
if parsedargs.ana_curv_h_opt:
    params["ana_is_optimally_curved_h"] = True
    params["ana_is_curved_h"] = True


if params["verbose"]:
    print("ki = %g / A, kf = %g / A, E = %g meV, Q = %g / A." %
        (params["ki"], params["kf"], params["E"], params["Q"]))
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# calculate resolution ellipsoid using the given backend
# -----------------------------------------------------------------------------
import reso
import pop
import eck

if params["reso_method"] == "eck":
    if params["verbose"]:
        print("\nCalculating using Eckold-Sobolev method. Scattering %s in kf." %
              ("vertically" if params["kf_vert"] else "horizontally"))
    res = eck.calc(params)
elif params["reso_method"] == "pop":
    if params["verbose"]:
        print("\nCalculating using Popovici method.")
    res = pop.calc(params, False)
elif params["reso_method"] == "cn":
    if params["verbose"]:
        print("\nsCalculating using Cooper-Nathans method.")
    res = pop.calc(params, True)
else:
    raise "ResPy: Invalid resolution calculation method selected."

if not res["ok"]:
    print("RESOLUTION CALCULATION FAILED!")
    exit(-1)


if params["verbose"]:
    print("R0 = %g, Vol = %g" % (res["r0"], res["res_vol"]))
    print("Resolution matrix:\n%s" % res["reso"])
    print("Resolution vector: %s" % res["reso_v"])
    print("Resolution scalar: %g" % res["reso_s"])
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# describe and plot ellipses
# -----------------------------------------------------------------------------
ellipses = reso.calc_ellipses(res["reso"], params["verbose"])


# plot ellipses
if show_plots:
    reso.plot_ellipses(ellipses, params["verbose"])
# -----------------------------------------------------------------------------
