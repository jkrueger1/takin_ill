/**
 * swig interface for magdyn
 * @author Tobias Weber <tweber@ill.fr>
 * @date 12-oct-2023
 * @license see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * mag-core (part of the Takin software suite)
 * Copyright (C) 2018-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
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

%module magdyn
%{
	#include "../../tlibs2/libs/magdyn.h"
%}


%include "std_vector.i"
%include "std_array.i"
%include "std_string.i"
%include "std_shared_ptr.i"
%include "std_complex.i"


%template(CplxD) std::complex<double>;

%template(VecD) std::vector<double>;
%template(VecCplx) std::vector<std::complex<double>>;

%template(ArrD3) std::array<double, 3>;
%template(ArrD4) std::array<double, 4>;
%template(ArrCplx3) std::array<std::complex<double>, 3>;
%template(ArrCplx4) std::array<std::complex<double>, 4>;
%template(ArrStr3) std::array<std::string, 3>;
%template(ArrStr33) std::array<std::array<std::string, 3>, 3>;


// ----------------------------------------------------------------------------
// matrix and vector containers
// ----------------------------------------------------------------------------
//%include "maths.h"

// TODO: add these as soon as swig supports concepts

//%template(VectorD) tl2::vec<double>;
//%template(VectorCplx) tl2::vec<std::complex<double>>;
//%template(MatrixD) tl2::mat<double>;
//%template(MatrixCplx) tl2::mat<std::complex<double>>;
// ----------------------------------------------------------------------------


%include "../../tlibs2/libs/magdyn.h"

// ----------------------------------------------------------------------------
// input- and output structs (and vectors of them)
// ----------------------------------------------------------------------------
%template(MagneticSite) tl2_mag::t_MagneticSite<
	tl2::mat<std::complex<double>>,
	tl2::vec<std::complex<double>>,
	tl2::vec<double>,
	std::size_t,
	double>;
%template(VecMagneticSite) std::vector<
	tl2_mag::t_MagneticSite<
		tl2::mat<std::complex<double>>,
		tl2::vec<std::complex<double>>,
		tl2::vec<double>,
		std::size_t,
		double>
	>;
%template(VecMagneticSitePtr) std::vector<
	const tl2_mag::t_MagneticSite<
		tl2::mat<std::complex<double>>,
		tl2::vec<std::complex<double>>,
		tl2::vec<double>,
		std::size_t,
		double>*
	>;

%template(ExchangeTerm) tl2_mag::t_ExchangeTerm<
	tl2::mat<std::complex<double>>,
	tl2::vec<std::complex<double>>,
	tl2::vec<double>,
	std::size_t,
	std::complex<double>>;
%template(VecExchangeTerm) std::vector<
	tl2_mag::t_ExchangeTerm<
		tl2::mat<std::complex<double>>,
		tl2::vec<std::complex<double>>,
		tl2::vec<double>,
		std::size_t,
		std::complex<double>>
	>;

%template(ExternalField) tl2_mag::t_ExternalField<
	tl2::vec<double>,
	double>;

%template(EnergyAndWeight) tl2_mag::t_EnergyAndWeight<
	tl2::mat<std::complex<double>>,
	double>;
%template(VecEnergyAndWeight) std::vector<
	tl2_mag::t_EnergyAndWeight<
		tl2::mat<std::complex<double>>,
		double>
	>;

%template(Variable) tl2_mag::t_Variable<
	std::complex<double>>;
%template(VecVariable) std::vector<
	tl2_mag::t_Variable<
		std::complex<double>>
	>;
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
/**
 * main magdyn class
 */
%template(MagDyn) tl2_mag::MagDyn<
	tl2::mat<std::complex<double>>,
	tl2::vec<std::complex<double>>,
	tl2::mat<double>,
	tl2::vec<double>,
	std::complex<double>,
	double,
	std::size_t>;
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
// additional functions
// ----------------------------------------------------------------------------
%inline
%{
	/**
	 * types
	 */
	using t_str = std::string;
	using t_real = double;
	using t_cplx = std::complex<double>;

	using t_MagDyn = tl2_mag::MagDyn<
		tl2::mat<std::complex<double>>,
		tl2::vec<std::complex<double>>,
		tl2::mat<double>, tl2::vec<double>,
		t_cplx, t_real, std::size_t>;


	/**
	 * adds a variable
	 */
	void add_variable(t_MagDyn& magdyn,
		const t_str& name,
		const t_cplx& val)
	{
		typename t_MagDyn::Variable var{};

		var.name = name;
		var.value = val;

		magdyn.AddVariable(std::move(var));
	}


	/**
	 * adds a magnetic site
	 */
	void add_site(t_MagDyn& magdyn,
		const t_str& name,
		const t_str& x, const t_str& y, const t_str& z,
		const t_str& sx = "0", const t_str& sy = "0", const t_str& sz = "1",
		const t_str& S = "1")
	{
		typename t_MagDyn::MagneticSite site{};

		site.name = name;

		site.pos[0] = x;
		site.pos[1] = y;
		site.pos[2] = z;

		site.spin_dir[0] = sx;
		site.spin_dir[1] = sy;
		site.spin_dir[2] = sz;
		site.spin_mag = S;

		magdyn.AddMagneticSite(std::move(site));
	}


	/**
	 * adds a coupling term between two magnetic sites
	 */
	void add_coupling(t_MagDyn& magdyn,
		const t_str& name,
		const t_str& site1, const t_str& site2,
		const t_str& dx, const t_str& dy, const t_str& dz,
		const t_str& J,
		const t_str& dmix = "", const t_str& dmiy = "", const t_str& dmiz = "",
		const t_str& Jxx = "", const t_str& Jxy = "", const t_str& Jxz = "",
		const t_str& Jyx = "", const t_str& Jyy = "", const t_str& Jyz = "",
		const t_str& Jzx = "", const t_str& Jzy = "", const t_str& Jzz = "")
	{
		typename t_MagDyn::ExchangeTerm coupling{};

		coupling.name = name;

		coupling.site1 = site1;
		coupling.site2 = site2;

		coupling.dist[0] = dx;
		coupling.dist[1] = dy;
		coupling.dist[2] = dz;

		coupling.J = J;

		coupling.dmi[0] = dmix;
		coupling.dmi[1] = dmiy;
		coupling.dmi[2] = dmiz;

		coupling.Jgen[0][0] = Jxx;
		coupling.Jgen[0][1] = Jxy;
		coupling.Jgen[0][2] = Jxz;
		coupling.Jgen[1][0] = Jyx;
		coupling.Jgen[1][1] = Jyy;
		coupling.Jgen[1][2] = Jyz;
		coupling.Jgen[2][0] = Jzx;
		coupling.Jgen[2][1] = Jzy;
		coupling.Jgen[2][2] = Jzz;

		magdyn.AddExchangeTerm(std::move(coupling));
	}


	/**
	 * calculates sites and couplings
	 */
	void calc(t_MagDyn& magdyn)
	{
		magdyn.CalcMagneticSites();
		magdyn.CalcExchangeTerms();
	}
%}
// ----------------------------------------------------------------------------
