#
# @author Tobias Weber <tweber@ill.fr>
# @date 2016
# @license GPLv2
#
# ----------------------------------------------------------------------------
# tlibs
# Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
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

find_library(Dl_LIBRARIES
	NAMES libdl.so libdl.so.2
	HINTS /usr/local/lib64 /usr/local/lib /usr/lib64 /usr/lib /opt/local/lib /lib/x86_64-linux-gnu #/usr/lib32 /usr/local/lib32
	DOC "DL library"
)

if(Dl_LIBRARIES)
	set(Dl_FOUND TRUE)
	message("DL library: ${Dl_LIBRARIES}")
else()
	set(Dl_FOUND FALSE)
	message("Error: libdl could not be found!")
endif()
