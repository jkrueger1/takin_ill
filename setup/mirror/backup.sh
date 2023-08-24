#!/bin/bash

#
# backup of source repositories
# @author Tobias Weber <tweber@ill.fr>
# @date mar-21
# @license GPLv2
#
# ----------------------------------------------------------------------------
# Takin (inelastic neutron scattering software package)
# Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
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

declare -a REPOS_ILL=(
	"https://github.com/ILLGrenoble/takin"
	"https://github.com/ILLGrenoble/taspaths"
	"https://code.ill.fr/scientific-software/takin/core"
	"https://code.ill.fr/scientific-software/takin/mag-core"
	"https://code.ill.fr/scientific-software/takin/tlibs"
	"https://code.ill.fr/scientific-software/takin/tlibs2"
	"https://code.ill.fr/scientific-software/takin/setup"
	"https://code.ill.fr/scientific-software/takin/data"
	"https://code.ill.fr/scientific-software/takin/plugins/mnsi"
	"https://code.ill.fr/scientific-software/takin/plugins/magnons"
	"https://code.ill.fr/scientific-software/takin/paths"
)


declare -a REPOS_FRM=(
	"https://forge.frm2.tum.de/cgit/cgit.cgi/frm2/mira/tastools.git"
	"https://forge.frm2.tum.de/cgit/cgit.cgi/frm2/mira/tlibs.git"
	"https://forge.frm2.tum.de/cgit/cgit.cgi/frm2/mira/takin-data.git"
	"https://forge.frm2.tum.de/cgit/cgit.cgi/frm2/mira/miezetools.git"
)


mkdir -v -p ill
pushd ill
for repo in ${REPOS_ILL[@]}; do
	echo "Cloning ${repo}..."
	git clone -v --recurse-submodules ${repo}
	echo -e "\n"
done
popd


mkdir -v -p frm
pushd frm
for repo in ${REPOS_FRM[@]}; do
	echo "Cloning ${repo}..."
	git clone -v --recurse-submodules ${repo}
	echo -e "\n"
done
popd
