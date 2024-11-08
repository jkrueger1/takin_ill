/**
 * calculates the berry curvatures
 * @author Tobias Weber <tweber@ill.fr>
 * @date 5-November-2024
 * @license see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
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

// g++ -std=c++20 -march=native -O2 -Wall -Wextra -Weffc++ -I .. -o berry berry.cpp -llapacke -llapack -lblas -lgfortran


#include "libs/magdyn.h"



// types
using t_real = double;
using t_cplx = std::complex<t_real>;
using t_mat = tl2::mat<t_cplx>;
using t_vec = tl2::vec<t_cplx>;
using t_mat_real = tl2::mat<t_real>;
using t_vec_real = tl2::vec<t_real>;
using t_magdyn = tl2_mag::MagDyn<
	t_mat, t_vec, t_mat_real, t_vec_real,
	t_cplx, t_real, std::size_t>;
using t_SofQE = typename t_magdyn::SofQE;


static constexpr t_real eps = 1e-4;
static constexpr unsigned int prec = 4;


using namespace tl2_ops;



void print_states(const t_SofQE& S)
{
	std::cout << "\nQ = " << S.Q_rlu << ", E = ";
	for(const auto& EandS : S.E_and_S)
		std::cout << EandS.E << ", ";

	std::cout << "states = \n";
	tl2::niceprint(std::cout, S.evec_mat, eps, prec);
	std::cout << std::endl;
}



int main(int argc, char** argv)
{
	std::cout.precision(prec);
	unsigned int width = prec * 3;
	t_real delta = 0.001;  // for differentiation
	t_real h = 0., k = 0., l = 0.;

	if(argc < 2)
	{
		std::cerr << "Please specify a magdyn file." << std::endl;
		return -1;
	}

	t_magdyn magdyn{};
	magdyn.SetEpsilon(eps);
	magdyn.SetUniteDegenerateEnergies(false);

	if(!magdyn.Load(argv[1]))
	{
		std::cerr << "Could not load model.";
		return -1;
	}

	std::cout << std::left << std::setw(width) << "# q" << " ";
	std::cout << std::left << std::setw(width) << "E_1" << " ";
	std::cout << std::left << std::setw(width) << "Re(b_1)" << " ";
	std::cout << std::left << std::setw(width) << "Im(b_1)" << " ";
	std::cout << std::left << std::setw(width) << "...";
	std::cout << std::endl;

	for(t_real q = 0.; q < 1.; q += 0.005)
	{
		std::cout << std::left << std::setw(width) << q << " ";
		t_vec_real Q = tl2::create<t_vec_real>({ h + q, k, l });

		t_SofQE S = magdyn.CalcEnergies(Q, false);
		//print_states(S);

		/*std::vector<t_vec> conns = magdyn.GetBerryConnections(Q, 0.001);
		for(const t_vec& conn : conns)
			std::cout << conn << std::endl;
		std::cout << std::endl;*/

		std::vector<t_cplx> curves = magdyn.GetBerryCurvatures(Q, delta);
		for(std::size_t band = 0; band < curves.size(); ++band)
		{
			const t_cplx& curve = curves[band];
			t_real E = S.E_and_S[band].E;
			std::cout << std::left << std::setw(width) << E << " ";

			std::cout << std::left << std::setw(width) << curve.real() << " ";
			std::cout << std::left << std::setw(width) << curve.imag() << " ";
		}
		std::cout << std::endl;
	}

	return 0;
}
