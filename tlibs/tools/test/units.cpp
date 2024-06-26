/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 *
 * ----------------------------------------------------------------------------
 * tlibs -- a physical-mathematical C++ template library
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) version 3.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ----------------------------------------------------------------------------
 */

#include <iostream>
#include "../phys/units.h"
#include "../phys/neutrons.h"
#include "../log/debug.h"

int main()
{
	std::cout.precision(8);

	tl::t_length_si<double> len = 1.*tl::t_meters<double>;
	std::cout << len << std::endl;

	tl::t_energy_si<double> E = 1.*tl::t_meV<double>;
	std::cout << E << std::endl;

	tl::t_energy_si<float> fE = 1.f*tl::t_meV<float>;
	std::cout << fE << std::endl;

	std::cout << "d: " << tl::t_KSQ2E<double> << std::endl;
	std::cout << "f: " << tl::t_KSQ2E<float> << std::endl;
	std::cout << "ld: " << tl::t_KSQ2E<long double> << std::endl;

	std::cout << tl::get_typename<decltype(tl::co::hbar)>() << std::endl;
	std::cout << tl::get_typename<tl::co::hbar_t>() << std::endl;

	std::cout << tl::get_m_n<double>() << std::endl;
	std::cout << tl::get_m_n<float>() << std::endl;

	std::cout << tl::get_kB<double>() << std::endl;
	std::cout << tl::get_kB<float>() << std::endl;


	bool bImag;
	std::cout << "e2k: " << tl::E2k(E, bImag) << std::endl;
	std::cout << "e2k, direct: " << tl::E2k_direct(E, bImag) << std::endl;
	std::cout << "float e2k: " << tl::E2k(fE, bImag) << std::endl;
	std::cout << "float e2k, direct: " << tl::E2k_direct(fE, bImag) << std::endl;


	tl::t_wavenumber_si<double> k = 1./tl::t_angstrom<double>;
	tl::t_wavenumber_si<float> fk = 1.f/tl::t_angstrom<float>;

	std::cout << "k2e: " << tl::k2E(k) << std::endl;
	std::cout << "k2e, direct: " << tl::k2E_direct(k) << std::endl;
	std::cout << "float k2e: " << tl::k2E(fk) << std::endl;
	std::cout << "float k2e, direct: " << tl::k2E_direct(fk) << std::endl;

	std::cout << "p: " << double(tl::get_r_e<double>() * (-tl::get_mu_n<double>()/tl::get_mu_N<double>()) * 0.5 /
		(tl::get_one_femtometer<double>())) << " fm" << std::endl;
	std::cout << "p: " << tl::get_mu_N<double>()*tl::get_mu_B<double>()*tl::get_gamma_n<double>() *
		2.*tl::get_m_n<double>()/(tl::get_hbar<double>()*tl::get_hbar<double>()) << std::endl;


	tl::t_energy_si<double> Ei = 5.*tl::t_meV<double>;
	tl::t_energy_si<double> dE = -1.*tl::t_meV<double>;
	tl::t_angle_si<double> tt = 90.*tl::t_degrees<double>;
	tl::t_wavenumber_si<double> Q = 2./tl::t_angstrom<double>;
	std::cout << tl::kinematic_plane(1, Ei, dE, tt) << std::endl;

	return 0;
}
