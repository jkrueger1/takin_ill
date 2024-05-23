#!/bin/bash
#
# takin brew build script
# @author Tobias Weber <tweber@ill.fr>
# @date sep-2020
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

# individual building steps
setup_buildenv=1
setup_externals=1
build_takin=1
build_takin2=1
build_plugins=1
build_py_modules=1
build_package=1

use_syspy=0


export MACOSX_DEPLOYMENT_TARGET=11.0

NUM_CORES=$(nproc)


# get root dir of takin repos
TAKIN_ROOT=$(dirname $0)/../..
cd "${TAKIN_ROOT}"
TAKIN_ROOT=$(pwd)
echo -e "Takin root dir: ${TAKIN_ROOT}"


if [ $setup_buildenv -ne 0 ]; then
	echo -e "\n================================================================================"
	echo -e "Setting up build environment..."
	echo -e "================================================================================\n"

	pushd "${TAKIN_ROOT}/setup"
		if ! ./build_mac/buildenv_brew.sh; then
			exit -1
		fi
	popd
fi


if [ $setup_externals -ne 0 ]; then
	echo -e "\n================================================================================"
	echo -e "Getting external dependencies..."
	echo -e "================================================================================\n"

	pushd "${TAKIN_ROOT}/core"
		rm -rf tmp
		if ! ../setup/externals/setup_externals.sh; then
			exit -1
		fi
		if ! ../setup/externals/get_3rdparty_licenses.sh; then
			exit -1
		fi
	popd

	pushd "${TAKIN_ROOT}/mag-core"
		if ! ../setup/externals/setup_externals_mag.sh; then
			exit -1
		fi
	popd
fi


if [ $build_takin -ne 0 ]; then
	echo -e "\n================================================================================"
	echo -e "Building main Takin binary..."
	echo -e "================================================================================\n"

	pushd "${TAKIN_ROOT}/core"
		../setup/build_general/clean.sh

		mkdir -p build
		cd build

		#if ! cmake -DDISABLE_INTERPROC_XSI=True -DDEBUG=False ..; then
		if ! cmake -DDEBUG=False ..; then
			echo -e "Failed configuring core package."
			exit -1
		fi

		if ! make -j${NUM_CORES}; then
			echo -e "Failed building core package."
			exit -1
		fi
	popd
fi


if [ $build_takin2 -ne 0 ]; then
	echo -e "\n================================================================================"
	echo -e "Building Takin 2 tools..."
	echo -e "================================================================================\n"

	pushd "${TAKIN_ROOT}/mag-core"
		rm -rf build
		mkdir -p build
		cd build

		if ! cmake -DCMAKE_BUILD_TYPE=Release -DONLY_BUILD_FINISHED=True ..; then
			echo -e "Failed configuring mag-core package."
			exit -1
		fi

		if ! make -j${NUM_CORES}; then
			echo -e "Failed building mag-core package."
			exit -1
		fi

		# copy tools to Takin main dir
		cp -v tools/cif2xml/takin_cif2xml "${TAKIN_ROOT}"/core/bin/
		cp -v tools/cif2xml/takin_findsg "${TAKIN_ROOT}"/core/bin/
		cp -v tools/pol/takin_pol "${TAKIN_ROOT}"/core/bin/
		cp -v tools/bz/takin_bz "${TAKIN_ROOT}"/core/bin/
		cp -v tools/structfact/takin_structfact "${TAKIN_ROOT}"/core/bin/
		cp -v tools/magstructfact/takin_magstructfact "${TAKIN_ROOT}"/core/bin/
		cp -v tools/magdyn/takin_magdyn "${TAKIN_ROOT}"/core/bin/
		cp -v tools/scanbrowser/takin_scanbrowser "${TAKIN_ROOT}"/core/bin/
		cp -v tools/magsgbrowser/takin_magsgbrowser "${TAKIN_ROOT}"/core/bin/
		cp -v tools/moldyn/takin_moldyn "${TAKIN_ROOT}"/core/bin/
	popd
fi


if [ $build_plugins -ne 0 ]; then
	echo -e "\n================================================================================"
	echo -e "Building Takin plugins..."
	echo -e "================================================================================\n"

	pushd "${TAKIN_ROOT}/magnon-plugin"
		rm -rf build
		mkdir -p build
		cd build

		if ! cmake -DCMAKE_BUILD_TYPE=Release ..; then
			echo -e "Failed configuring magnon plugin."
			exit -1
		fi

		if ! make -j${NUM_CORES}; then
			echo -e "Failed building magnon plugin."
			exit -1
		fi

		# copy plugin to Takin main dir
		cp -v libmagnonmod.dylib "${TAKIN_ROOT}"/core/plugins/
	popd
fi


if [ $build_py_modules -ne 0 ]; then
	echo -e "\n================================================================================"
	echo -e "Building py modules..."
	echo -e "================================================================================\n"

	pushd "${TAKIN_ROOT}/mag-core/tools_py/magdyn"
		rm -rf build
		mkdir -p build
		cd build

		if ! cmake -DCMAKE_BUILD_TYPE=Release ..; then
			echo -e "Failed configuring magnetic dynamics py module."
			exit -1
		fi

		if ! make -j${NUM_CORES}; then
			echo -e "Failed building magnetic dynamics py module."
			exit -1
		fi

		cp -v _magdyn_py.dylib "${TAKIN_ROOT}"/core/pymods/
		cp -v _magdyn_py.so "${TAKIN_ROOT}"/core/pymods/
		cp -v magdyn.py "${TAKIN_ROOT}"/core/pymods/
	popd

	pushd "${TAKIN_ROOT}/mag-core/tools_py/instr"
		rm -rf build
		mkdir -p build
		cd build

		if ! cmake -DCMAKE_BUILD_TYPE=Release ..; then
			echo -e "Failed configuring instrument data loader py module."
			exit -1
		fi

		if ! make -j${NUM_CORES}; then
			echo -e "Failed building instrument data loader py module."
			exit -1
		fi

		cp -v _instr_py.dylib "${TAKIN_ROOT}"/core/pymods/
		cp -v _instr_py.so "${TAKIN_ROOT}"/core/pymods/
		cp -v instr.py "${TAKIN_ROOT}"/core/pymods/
	popd
fi



if [ $build_package -ne 0 ]; then
	pushd "${TAKIN_ROOT}"
		rm -rf tmp
		mkdir -p tmp

		cd core

		echo -e "\n================================================================================"
		echo -e "Creating icons..."
		echo -e "================================================================================\n"
		if ! ../setup/build_mac/mk_icon.sh; then
			exit -1
		fi

		echo -e "\n================================================================================"
		echo -e "Creating app directory..."
		echo -e "================================================================================\n"
		if ! ../setup/build_mac/cp_app.sh; then
			exit -1
		fi

		echo -e "\n================================================================================"
		echo -e "Cleaning up app directory..."
		echo -e "================================================================================\n"
		if ! ../setup/build_mac/clean_app.sh; then
			exit -1
		fi

		echo -e "\n================================================================================"
		echo -e "Fixing dynamic binding for local libraries..."
		echo -e "================================================================================\n"
		if ! ../setup/build_mac/fix_names.sh; then
			exit -1
		fi

		if [ $use_syspy -ne 0 ]; then
			echo -e "\n================================================================================"
			echo -e "Using system python frameworks instead..."
			echo -e "================================================================================\n"
			if ! ../setup/build_mac/use_syspy.sh; then
				exit -1
			fi
			if ! ../setup/build_mac/clean_app.sh; then
				exit -1
			fi
		fi

		echo -e "\n================================================================================"
		echo -e "Building Takin package..."
		echo -e "================================================================================\n"
		if ! ../setup/build_mac/cp_dmg.sh; then
			exit -1
		fi

		cp -v takin.dmg ../tmp/takin.dmg
	popd


	if [ -e  "${TAKIN_ROOT}/tmp/takin.dmg" ]; then
		echo -e "\n================================================================================"
		echo -e "The built Takin package can be found here:\n\t${TAKIN_ROOT}/tmp/takin.dmg"
		echo -e "================================================================================\n"
	else
		echo -e "\n================================================================================"
		echo -e "Error: Takin package could not be built!"
		echo -e "================================================================================\n"
	fi
fi
