/**
 * tlibs2 -- magnetic dynamics
 * @author Tobias Weber <tweber@ill.fr>
 * @date 2022 - 2024
 * @license GPLv3, see 'LICENSE' file
 *
 * References:
 *   - (Toth 2015) S. Toth and B. Lake, J. Phys.: Condens. Matter 27 166002 (2015):
 *                 https://doi.org/10.1088/0953-8984/27/16/166002
 *                 https://arxiv.org/abs/1402.6069
 *   - (Heinsdorf 2021) N. Heinsdorf, manual example calculation for a simple
 *                      ferromagnetic case, personal communications, 2021/2022.
 *
 * @desc This file implements the formalism given by (Toth 2015).
 * @desc For further references, see the 'LITERATURE' file.
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

#ifndef __TLIBS2_MAGDYN_H__
#define __TLIBS2_MAGDYN_H__

#include <algorithm>
#include <numeric>
#include <vector>
#include <array>
#include <tuple>
#include <unordered_set>
#include <unordered_map>
#include <string>
#include <string_view>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdint>

#include <boost/asio.hpp>
#include <boost/container_hash/hash.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#ifndef USE_LAPACK
	#define USE_LAPACK 1
#endif
#include "maths.h"
#include "units.h"
#include "phys.h"
#include "algos.h"
#include "expr.h"
#include "fit.h"

// enables debug output
//#define __TLIBS2_MAGDYN_DEBUG_OUTPUT__

// enables ground state minimisation
//#define __TLIBS2_MAGDYN_USE_MINUIT__



namespace tl2_mag {

// ----------------------------------------------------------------------------
// helper functions
// ----------------------------------------------------------------------------

/**
 * rotate spin vector for incommensurate structures, i.e. helices
 */
template<class t_mat, class t_vec, class t_real = typename t_mat::value_type>
#ifndef SWIG  // TODO: remove this as soon as swig understands concepts
requires tl2::is_mat<t_mat> && tl2::is_vec<t_vec>
#endif
void rotate_spin_incommensurate(t_vec& spin_vec,
	const t_vec& sc_vec, const t_vec& ordering, const t_vec& rotaxis,
	t_real eps = std::numeric_limits<t_real>::epsilon())
{
	const t_real sc_angle = t_real(2) * tl2::pi<t_real> *
		tl2::inner<t_vec>(ordering, sc_vec);

	if(!tl2::equals_0<t_real>(sc_angle, t_real(eps)))
	{
		t_mat sc_rot = tl2::rotation<t_mat, t_vec>(rotaxis, sc_angle);
		spin_vec = sc_rot * spin_vec;
	}
}



/**
 * create a 3-vector from a homogeneous 4-vector
 */
template<class t_vec>
#ifndef SWIG  // TODO: remove this as soon as swig understands concepts
requires tl2::is_vec<t_vec>
#endif
t_vec to_3vec(const t_vec& vec)
{
	return tl2::create<t_vec>({ vec[0], vec[1], vec[2] });
}



/**
 * create a (homogeneous) 4-vector from a 3-vector
 */
template<class t_vec, class t_val = typename t_vec::value_type>
#ifndef SWIG  // TODO: remove this as soon as swig understands concepts
requires tl2::is_vec<t_vec>
#endif
t_vec to_4vec(const t_vec& vec, t_val w = 0.)
{
	return tl2::create<t_vec>({ vec[0], vec[1], vec[2], w });
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// input- and output struct templates
// ----------------------------------------------------------------------------

using t_strarr3 = std::array<std::string, 3>;
using t_strarr33 = std::array<std::array<std::string, 3>, 3>;



/**
 * magnetic sites
 */
template<class t_mat, class t_vec, class t_vec_real,
	class t_size, class t_real = typename t_vec_real::value_type>
#ifndef SWIG  // TODO: remove this as soon as swig understands concepts
requires tl2::is_mat<t_mat> && tl2::is_vec<t_vec> && tl2::is_vec<t_vec_real>
#endif
struct t_MagneticSite
{
	// ------------------------------------------------------------------------
	// input properties
	std::string name{};          // identifier
	t_size sym_idx{};            // groups positions belonging to the same symmetry group (0: none)

	t_strarr3 pos{};             // magnetic site position

	t_strarr3 spin_dir{};        // spin direction
	t_strarr3 spin_ortho{};      // spin orthogonal vector

	std::string spin_mag{};      // spin magnitude
	t_mat g_e{};                 // electron g factor
	// ------------------------------------------------------------------------

	// ------------------------------------------------------------------------
	// calculated properties
	t_vec_real pos_calc{};       // magnetic site position

	t_vec_real spin_dir_calc{};  // spin vector
	t_vec trafo_z_calc{};        // trafo z vector (3rd column in trafo matrix)
	t_vec trafo_plane_calc{};    // trafo orthogonal plane (1st and 2nd coumns)
	t_vec trafo_plane_conj_calc{};

	t_vec ge_trafo_z_calc{};     // g_e * trafo z vector
	t_vec ge_trafo_plane_calc{}; // g_e * trafo orthogonal plane
	t_vec ge_trafo_plane_conj_calc{};

	t_real spin_mag_calc{};      // spin magnitude
	// ------------------------------------------------------------------------
};



/**
 * couplings between magnetic sites
 */
template<class t_mat, class t_vec, class t_vec_real,
	class t_size, class t_cplx = typename t_mat::value_type,
	class t_real = typename t_vec_real::value_type>
#ifndef SWIG  // TODO: remove this as soon as swig understands concepts
requires tl2::is_mat<t_mat> && tl2::is_vec<t_vec>
#endif
struct t_ExchangeTerm
{
	// ------------------------------------------------------------------------
	// input properties (parsable expressions)
	std::string name{};          // identifier
	t_size sym_idx{};            // groups couplings belonging to the same symmetry group (0: none)

	std::string site1{};         // name of first magnetic site
	std::string site2{};         // name of second magnetic site
	t_strarr3 dist{};            // distance between unit cells

	std::string J{};             // Heisenberg interaction
	t_strarr3 dmi{};             // Dzyaloshinskij-Moriya interaction
	t_strarr33 Jgen{};           // general exchange interaction
	// ------------------------------------------------------------------------

	// ------------------------------------------------------------------------
	// calculated properties
	t_size site1_calc{};         // index of first magnetic site
	t_size site2_calc{};         // index of second magnetic site
	t_vec_real dist_calc{};      // distance between unit cells
	t_real length_calc{};        // length of the coupling (in lab units)

	t_cplx J_calc{};             // Heisenberg interaction
	t_vec dmi_calc{};            // Dzyaloshinskij-Moriya interaction
	t_mat Jgen_calc{};           // general exchange interaction
	// ------------------------------------------------------------------------
};



/**
 * terms related to an external magnetic field
 */
template<class t_vec_real, class t_real>
#ifndef SWIG  // TODO: remove this as soon as swig understands concepts
requires tl2::is_vec<t_vec_real>
#endif
struct t_ExternalField
{
	bool align_spins{};          // align spins along external field

	t_vec_real dir{};            // field direction
	t_real mag{};                // field magnitude
};



/**
 * eigenenergies and spin-spin correlation matrix
 */
template<class t_mat, class t_real, class t_cplx = typename t_mat::value_type>
#ifndef SWIG  // TODO: remove this as soon as swig understands concepts
requires tl2::is_mat<t_mat>
#endif
struct t_EnergyAndWeight
{
	t_real E{};

	// full dynamical structure factor
	t_mat S{};
	t_cplx S_sum{};
	t_real weight_full{};

	// projected dynamical structure factor for neutron scattering
	t_mat S_perp{};
	t_cplx S_perp_sum{};
	t_real weight{};
};



template<class t_mat, class t_real, class t_cplx = typename t_mat::value_type>
#ifndef SWIG  // TODO: remove this as soon as swig understands concepts
requires tl2::is_mat<t_mat>
#endif
struct t_SofQE
{
	t_real h, k, l;
	std::vector<t_EnergyAndWeight<t_mat, t_real, t_cplx>> E_and_S;
};



/**
 * variables for the expression parser
 */
template<class t_cplx>
struct t_Variable
{
	std::string name{};
	t_cplx value{};
};
// ----------------------------------------------------------------------------



/**
 * calculates magnon dynamics,
 * implementing the formalism given in (Toth 2015)
 */
template<
	class t_mat, class t_vec,
	class t_mat_real, class t_vec_real,
	class t_cplx = typename t_mat::value_type,
	class t_real = typename t_mat_real::value_type,
	class t_size = std::size_t>
#ifndef SWIG  // TODO: remove this as soon as swig understands concepts
requires tl2::is_mat<t_mat> && tl2::is_vec<t_vec> &&
	tl2::is_mat<t_mat_real> && tl2::is_vec<t_vec_real>
#endif
class MagDyn
{
public:
	// --------------------------------------------------------------------
	// structs and types
	// --------------------------------------------------------------------
	using MagneticSite = t_MagneticSite<t_mat, t_vec, t_vec_real, t_size, t_real>;
	using MagneticSites = std::vector<MagneticSite>;

	using ExchangeTerm = t_ExchangeTerm<t_mat, t_vec, t_vec_real, t_size, t_cplx, t_real>;
	using ExchangeTerms = std::vector<ExchangeTerm>;

	using Variable = t_Variable<t_cplx>;
	using Variables = std::vector<Variable>;

	using ExternalField = t_ExternalField<t_vec_real, t_real>;

	using EnergyAndWeight = t_EnergyAndWeight<t_mat, t_real, t_cplx>;
	using EnergiesAndWeights = std::vector<EnergyAndWeight>;

	using SofQE = t_SofQE<t_mat, t_real, t_cplx>;
	using SofQEs = std::vector<SofQE>;

	using t_indices = std::pair<t_size, t_size>;
	using t_Jmap = std::unordered_map<t_indices, t_mat, boost::hash<t_indices>>;
	// --------------------------------------------------------------------



public:
	MagDyn() = default;
	~MagDyn() = default;

	MagDyn(const MagDyn&) = default;
	MagDyn& operator=(const MagDyn&) = default;



	// --------------------------------------------------------------------
	// cleanup functions
	// --------------------------------------------------------------------
	/**
	 * clear all
	 */
	void Clear()
	{
		ClearVariables();
		ClearMagneticSites();
		ClearExchangeTerms();
		ClearExternalField();

		// clear temperature, -1: don't use
		m_temperature = -1.;

		// clear form factor
		m_magffact_formula = "";

		// clear ordering wave vector
		m_ordering = tl2::zero<t_vec_real>(3);
		m_is_incommensurate = false;

		// reset rotation axis
		m_rotaxis = tl2::create<t_vec_real>({ 1., 0., 0. });

		// clear crystal
		m_xtallattice[0] = m_xtallattice[1] = m_xtallattice[2] = 0.;
		m_xtalangles[0] = m_xtalangles[1] = m_xtalangles[2] = t_real(0.5) * tl2::pi<t_real>;
		m_xtalA = m_xtalB = tl2::unit<t_mat_real>(3);
		m_xtalUB = m_xtalUBinv = tl2::unit<t_mat_real>(3);

		// clear scattering plane
		m_scatteringplane[0] = tl2::create<t_vec_real>({ 1., 0., 0. });
		m_scatteringplane[1] = tl2::create<t_vec_real>({ 0., 1., 0. });
		m_scatteringplane[2] = tl2::create<t_vec_real>({ 0., 0., 1. });
	}



	/**
	 * clear all parser variables
	 */
	void ClearVariables()
	{
		m_variables.clear();
	}



	/**
	 * clear all magnetic sites
	 */
	void ClearMagneticSites()
	{
		m_sites.clear();
	}



	/**
	 * clear all couplings
	 */
	void ClearExchangeTerms()
	{
		m_exchange_terms.clear();
	}



	/**
	 * clear the external field settings
	 */
	void ClearExternalField()
	{
		m_field.dir.clear();
		m_field.mag = 0.;
		m_field.align_spins = false;
	}
	// --------------------------------------------------------------------



	// --------------------------------------------------------------------
	// getter
	// --------------------------------------------------------------------
	const Variables& GetVariables() const { return m_variables; }
	const MagneticSites& GetMagneticSites() const { return m_sites; }
	MagneticSites& GetMagneticSites() { return m_sites; }
	t_size GetMagneticSitesCount() const { return m_sites.size(); }
	const ExchangeTerms& GetExchangeTerms() const { return m_exchange_terms; }
	ExchangeTerms& GetExchangeTerms() { return m_exchange_terms; }
	t_size GetExchangeTermsCount() const { return m_exchange_terms.size(); }

	const ExternalField& GetExternalField() const { return m_field; }
	const t_vec_real& GetRotationAxis() const { return m_rotaxis; }
	const t_vec_real& GetOrderingWavevector() const { return m_ordering; }

	t_real GetTemperature() const { return m_temperature; }
	t_real GetBoseCutoffEnergy() const { return m_bose_cutoff; }

	const std::string& GetMagneticFormFactor() const { return m_magffact_formula; }



	const MagneticSite& GetMagneticSite(t_size idx) const
	{
		if(!CheckMagneticSite(idx))
		{
			static MagneticSite null_site{};
			return null_site;
		}

		return m_sites[idx];
	}



	const ExchangeTerm& GetExchangeTerm(t_size idx) const
	{
		if(!CheckExchangeTerm(idx))
		{
			static ExchangeTerm null_term{};
			return null_term;
		}

		return m_exchange_terms[idx];
	}



	bool IsIncommensurate() const
	{
		return m_is_incommensurate || m_force_incommensurate;
	}



	/**
	 * get number of magnetic sites with the given name (to check if the name is unique)
	 */
	std::vector<const MagneticSite*> FindMagneticSites(const std::string& name) const
	{
		std::vector<const MagneticSite*> sites;

		for(const MagneticSite& site : GetMagneticSites())
		{
			if(site.name == name)
				sites.push_back(&site);
		}

		return sites;
	}



	/**
	 * get magnetic site with the given name
	 */
	const MagneticSite* FindMagneticSite(const std::string& name) const
	{
		for(const MagneticSite& site : GetMagneticSites())
		{
			if(site.name == name)
				return &site;
		}

		return nullptr;
	}



	/**
	 * get the index of a magnetic site from its name
	 */
	t_size GetMagneticSiteIndex(const std::string& name) const
	{
		// try to find the site index by name
		for(t_size idx = 0; idx < GetMagneticSitesCount(); ++idx)
		{
			if(GetMagneticSite(idx).name == name)
				return idx;
		}

		// alternatively try to parse the expression for the index
		tl2::ExprParser<t_size> parser;
		parser.SetInvalid0(false);
		parser.SetAutoregisterVariables(false);
		if(parser.parse_noexcept(name))
		{
			if(t_size idx = parser.eval_noexcept(); idx < GetMagneticSitesCount())
				return idx;
		}
		else
		{
			std::cerr << "Magdyn error: "
				<< "Invalid site name \"" << name << "\"."
				<< std::endl;
		}

		// nothing found: return invalid index
		return GetMagneticSitesCount();
	}



	/**
	 * get the index of an exchange term from its name
	 */
	t_size GetExchangeTermIndex(const std::string& name) const
	{
		// try to find the term index by name
		for(t_size idx = 0; idx < GetExchangeTermsCount(); ++idx)
		{
			if(GetExchangeTerm(idx).name == name)
				return idx;
		}

		// alternatively try to parse the expression for the index
		tl2::ExprParser<t_size> parser;
		parser.SetInvalid0(false);
		parser.SetAutoregisterVariables(false);
		if(parser.parse_noexcept(name))
		{
			if(t_size idx = parser.eval_noexcept(); idx < GetExchangeTermsCount())
				return idx;
		}
		else
		{
			std::cerr << "Magdyn error: Invalid coupling name \"" << name << "\"."
				<< std::endl;
		}

		// nothing found: return invalid index
		return GetExchangeTermsCount();
	}



	std::vector<t_vec_real> GetMagneticSitePositions(bool homogeneous = false) const
	{
		std::vector<t_vec_real> sites;
		sites.reserve(GetMagneticSitesCount());

		for(const MagneticSite& site : GetMagneticSites())
		{
			if(homogeneous)
				sites.push_back(to_4vec<t_vec_real>(site.pos_calc, 1.));
			else
				sites.push_back(to_3vec<t_vec_real>(site.pos_calc));
		}

		return sites;
	}



	t_vec_real GetCrystalLattice() const
	{
		return tl2::create<t_vec_real>(
		{
			m_xtallattice[0], m_xtallattice[1], m_xtallattice[2],
			m_xtalangles[0], m_xtalangles[1], m_xtalangles[2],
		});
	}



	const t_vec_real* GetScatteringPlane() const
	{
		return m_scatteringplane;
	}



	/**
	 * get the needed supercell ranges from the exchange terms
	 */
	std::tuple<t_vec_real, t_vec_real> GetSupercellMinMax() const
	{
		t_vec_real min = tl2::zero<t_vec_real>(3);
		t_vec_real max = tl2::zero<t_vec_real>(3);

		for(const ExchangeTerm& term : GetExchangeTerms())
		{
			for(std::uint8_t i = 0; i < 3; ++i)
			{
				min[i] = std::min(min[i], term.dist_calc[i]);
				max[i] = std::max(max[i], term.dist_calc[i]);
			}
		}

		return std::make_tuple(min, max);
	}
	// --------------------------------------------------------------------



	// --------------------------------------------------------------------
	// setter
	// --------------------------------------------------------------------
	void SetEpsilon(t_real eps) { m_eps = eps; }
	void SetPrecision(int prec) { m_prec = prec; }

	void SetTemperature(t_real T) { m_temperature = T; }
	void SetBoseCutoffEnergy(t_real E) { m_bose_cutoff = E; }
	void SetUniteDegenerateEnergies(bool b) { m_unite_degenerate_energies = b; }
	void SetForceIncommensurate(bool b) { m_force_incommensurate = b; }
	void SetPerformChecks(bool b) { m_perform_checks = b; }

	void SetPhaseSign(t_real sign) { m_phase_sign = sign; }
	void SetCholeskyMaxTries(t_size max_tries) { m_tries_chol = max_tries; }
	void SetCholeskyInc(t_real delta) { m_delta_chol = delta; }



	void SetMagneticFormFactor(const std::string& ffact)
	{
		m_magffact_formula = ffact;
		if(m_magffact_formula == "")
			return;

		// parse the given formula
		m_magffact = GetExprParser();
		m_magffact.SetInvalid0(false);
		m_magffact.register_var("Q", 0.);

		if(!m_magffact.parse_noexcept(ffact))
		{
			m_magffact_formula = "";

			std::cerr << "Magdyn error: Magnetic form facor formula: \""
				<< ffact << "\" could not be parsed."
				<< std::endl;
		}
	}



	void SetExternalField(const ExternalField& field)
	{
		m_field = field;

		// normalise direction vector
		const t_real len = tl2::norm<t_vec_real>(m_field.dir);
		if(!tl2::equals_0<t_real>(len, m_eps))
			m_field.dir /= len;
	}



	void RotateExternalField(const t_vec_real& axis, t_real angle)
	{
		const t_mat_real rot = tl2::rotation<t_mat_real, t_vec_real>(
			axis, angle, false);
		m_field.dir = rot * m_field.dir;
	}



	void RotateExternalField(t_real x, t_real y, t_real z, t_real angle)
	{
		RotateExternalField(tl2::create<t_vec_real>({ x, y, z }), angle);
	}



	/**
	 * set the ordering wave vector (e.g., the helix pitch) for incommensurate structures
	 */
	void SetOrderingWavevector(const t_vec_real& ordering)
	{
		m_ordering = ordering;
		m_is_incommensurate = !tl2::equals_0<t_vec_real>(m_ordering, m_eps);
	}



	/**
	 * set the rotation axis for the ordering wave vector
	 */
	void SetRotationAxis(const t_vec_real& axis)
	{
		m_rotaxis = axis;

		// normalise
		const t_real len = tl2::norm<t_vec_real>(m_rotaxis);
		if(!tl2::equals_0<t_real>(len, m_eps))
			m_rotaxis /= len;
	}



	void SetCalcHamiltonian(bool H, bool Hp, bool Hm)
	{
		m_calc_H  = H;
		m_calc_Hp = Hp;
		m_calc_Hm = Hm;
	}



	void AddVariable(Variable&& var)
	{
		m_variables.emplace_back(std::forward<Variable&&>(var));
	}



	void SetVariable(Variable&& var)
	{
		// is a variable with the same name already registered?
		auto iter = std::find_if(m_variables.begin(), m_variables.end(),
			[&var](const Variable& thevar)
		{
			return thevar.name == var.name;
		});

		if(iter == m_variables.end())
		{
			// add a new variable
			AddVariable(std::forward<Variable&&>(var));
		}
		else
		{
			// replace the value of an existing variable
			iter->value = var.value;
		}
	}



	void AddMagneticSite(MagneticSite&& site)
	{
		m_sites.emplace_back(std::forward<MagneticSite&&>(site));
	}



	void AddExchangeTerm(ExchangeTerm&& term)
	{
		m_exchange_terms.emplace_back(std::forward<ExchangeTerm&&>(term));
	}



	/**
	 * calculate the B matrix from the crystal lattice
	 */
	void SetCrystalLattice(t_real a, t_real b, t_real c,
		t_real alpha, t_real beta, t_real gamma)
	{
		try
		{
			m_xtallattice[0] = a;
			m_xtallattice[1] = b;
			m_xtallattice[2] = c;

			m_xtalangles[0] = alpha;
			m_xtalangles[1] = beta;
			m_xtalangles[2] = gamma;

			// crystal fractional coordinate trafo matrices
			m_xtalA = tl2::A_matrix<t_mat_real>(a, b, c, alpha, beta, gamma);
			m_xtalB = tl2::B_matrix<t_mat_real>(a, b, c, alpha, beta, gamma);
		}
		catch(const std::exception& ex)
		{
			m_xtalA = m_xtalB = tl2::unit<t_mat_real>(3);
			std::cerr << "Magdyn error: Could not calculate crystal matrices."
				<< std::endl;
		}
	}



	/**
	 * calculate the UB matrix from the scattering plane and the crystal lattice
	 * note: SetCrystalLattice() has to be called before this function
	 */
	void SetScatteringPlane(t_real ah, t_real ak, t_real al,
		t_real bh, t_real bk, t_real bl)
	{
		try
		{
			m_scatteringplane[0] = tl2::create<t_vec_real>({ ah, ak, al });
			m_scatteringplane[1] = tl2::create<t_vec_real>({ bh, bk, bl });
			m_scatteringplane[2] = tl2::cross(m_xtalB, m_scatteringplane[0], m_scatteringplane[1]);

			m_xtalUB = tl2::UB_matrix(m_xtalB, m_scatteringplane[0], m_scatteringplane[1], m_scatteringplane[2]);
			bool inv_ok = false;
			std::tie(m_xtalUBinv, inv_ok) = tl2::inv(m_xtalUB);

			if(!inv_ok)
				std::cerr << "Magdyn error: UB matrix is not invertible."
					<< std::endl;
		}
		catch(const std::exception& ex)
		{
			m_xtalUB = m_xtalUBinv = tl2::unit<t_mat_real>(3);
			std::cerr << "Magdyn error: Could not calculate scattering plane matrices."
				<< std::endl;
		}
	}
	// --------------------------------------------------------------------



	// --------------------------------------------------------------------
	/**
	 * get an expression parser object with registered variables
	 */
	tl2::ExprParser<t_cplx> GetExprParser() const
	{
		tl2::ExprParser<t_cplx> parser;

		// register all variables
		parser.SetAutoregisterVariables(false);
		for(const Variable& var : m_variables)
			parser.register_var(var.name, var.value);

		return parser;
	}
	// --------------------------------------------------------------------



	// --------------------------------------------------------------------
	// sanity checks
	// --------------------------------------------------------------------
	/**
	 * check if the site index is valid
	 */
	bool CheckMagneticSite(t_size idx, bool print_error = true) const
	{
		if(!m_perform_checks)
			return true;

		if(idx >= m_sites.size())
		{
			if(print_error)
			{
				std::cerr << "Magdyn error: Site index " << idx
					<< " is out of bounds."
					<< std::endl;
			}

			return false;
		}

		return true;
	}



	/**
	 * check if the term index is valid
	 */
	bool CheckExchangeTerm(t_size idx, bool print_error = true) const
	{
		if(!m_perform_checks)
			return true;

		if(idx >= m_exchange_terms.size())
		{
			if(print_error)
			{
				std::cerr << "Magdyn error: Coupling index " << idx
					<< " is out of bounds."
					<< std::endl;
			}

			return false;
		}

		return true;
	}



	/**
	 * check if imaginary weights remain
	 */
	bool CheckImagWeights(const t_vec_real& Q_rlu, const EnergiesAndWeights& Es_and_S) const
	{
		if(!m_perform_checks)
			return true;

		using namespace tl2_ops;
		bool ok = true;

		for(const EnergyAndWeight& EandS : Es_and_S)
		{
			// imaginary parts should be gone after UniteEnergies()
			if(std::abs(EandS.S_perp_sum.imag()) > m_eps ||
				std::abs(EandS.S_sum.imag()) > m_eps)
			{
				ok = false;

				std::cerr << "Magdyn warning: Remaining imaginary S(Q, E) component at Q = "
					<< Q_rlu << " and E = " << EandS.E
					<< ": imag(S) = " << EandS.S_sum.imag()
					<< ", imag(S_perp) = " << EandS.S_perp_sum.imag()
					<< "." << std::endl;
			}
		}

		return ok;
	}
	// --------------------------------------------------------------------



	// --------------------------------------------------------------------
	// symmetrisation and generation functions
	// --------------------------------------------------------------------
	/**
	 * generate symmetric positions based on the given symops
	 */
	void SymmetriseMagneticSites(const std::vector<t_mat_real>& symops)
	{
		CalcExternalField();
		CalcMagneticSites();

		MagneticSites newsites;
		newsites.reserve(GetMagneticSitesCount() * symops.size());

		for(const MagneticSite& site : GetMagneticSites())
		{
			// get symmetry-equivalent positions
			const auto positions = tl2::apply_ops_hom<t_vec_real, t_mat_real, t_real>(
				site.pos_calc, symops, m_eps);

			for(t_size idx = 0; idx < positions.size(); ++idx)
			{
				MagneticSite newsite = site;
				newsite.pos_calc = positions[idx];
				newsite.pos[0] = tl2::var_to_str(newsite.pos_calc[0], m_prec);
				newsite.pos[1] = tl2::var_to_str(newsite.pos_calc[1], m_prec);
				newsite.pos[2] = tl2::var_to_str(newsite.pos_calc[2], m_prec);
				newsite.name += "_" + tl2::var_to_str(idx + 1, m_prec);

				newsites.emplace_back(std::move(newsite));
			}
		}

		m_sites = std::move(newsites);
		RemoveDuplicateMagneticSites();
		CalcSymmetryIndices(symops);
		CalcMagneticSites();
	}



	/**
	 * generate symmetric exchange terms based on the given symops
	 */
	void SymmetriseExchangeTerms(const std::vector<t_mat_real>& symops)
	{
		CalcExternalField();
		CalcMagneticSites();
		CalcExchangeTerms();

		ExchangeTerms newterms;
		newterms.reserve(GetExchangeTermsCount() * symops.size());

		tl2::ExprParser<t_cplx> parser = GetExprParser();

		// create unit cell site vectors
		std::vector<t_vec_real> sites_uc = GetMagneticSitePositions(true);

		for(const ExchangeTerm& term : GetExchangeTerms())
		{
			// check if the site indices are valid
			if(!CheckMagneticSite(term.site1_calc) || !CheckMagneticSite(term.site2_calc))
				continue;

			// super cell distance vector
			t_vec_real dist_sc = to_4vec<t_vec_real>(term.dist_calc, 0.);

			// generate new (possibly supercell) sites with symop
			auto sites1_sc = tl2::apply_ops_hom<t_vec_real, t_mat_real, t_real>(
				sites_uc[term.site1_calc], symops, m_eps,
				false /*keep in uc*/, true /*ignore occupied*/,
				true /*return homogeneous*/);
			auto sites2_sc = tl2::apply_ops_hom<t_vec_real, t_mat_real, t_real>(
				sites_uc[term.site2_calc] + dist_sc, symops, m_eps,
				false /*keep in uc*/, true /*ignore occupied*/,
				true /*return homogeneous*/);

			// generate new dmi vectors
			t_vec_real dmi = tl2::zero<t_vec_real>(4);

			for(std::uint8_t dmi_idx = 0; dmi_idx < 3; ++dmi_idx)
			{
				if(term.dmi[dmi_idx] == "")
					continue;

				if(parser.parse_noexcept(term.dmi[dmi_idx]))
				{
					dmi[dmi_idx] = parser.eval_noexcept().real();
				}
				else
				{
					std::cerr << "Magdyn error: Parsing DMI component " << dmi_idx
						<< " of term \"" << term.name << "\"."
						<< std::endl;
				}
			}

			const auto newdmis = tl2::apply_ops_hom<t_vec_real, t_mat_real, t_real>(
				dmi, symops, m_eps,
				false /*keep in uc*/, true /*ignore occupied*/,
				false /*return homogeneous*/, true /*pseudovector*/);

			// generate new general J matrices
			t_real Jgen_arr[3][3]{};

			for(std::uint8_t J_idx1 = 0; J_idx1 < 3; ++J_idx1)
			{
				for(std::uint8_t J_idx2 = 0; J_idx2 < 3; ++J_idx2)
				{
					if(term.Jgen[J_idx1][J_idx2] == "")
						continue;

					if(parser.parse_noexcept(term.Jgen[J_idx1][J_idx2]))
					{
						Jgen_arr[J_idx1][J_idx2] = parser.eval_noexcept().real();
					}
					else
					{
						std::cerr << "Magdyn error: Parsing general J component ("
							<< J_idx1 << ", " << J_idx2
							<< ") of term \"" << term.name << "\"."
							<< std::endl;
					}
				}
			}

			t_mat_real Jgen = tl2::create<t_mat_real>({
				Jgen_arr[0][0], Jgen_arr[0][1], Jgen_arr[0][2], 0,
				Jgen_arr[1][0], Jgen_arr[1][1], Jgen_arr[1][2], 0,
				Jgen_arr[2][0], Jgen_arr[2][1], Jgen_arr[2][2], 0,
				0,              0,              0,              0 });

			const auto newJgens = tl2::apply_ops_hom<t_mat_real, t_real>(Jgen, symops);

			// iterate and insert generated couplings
			for(t_size op_idx = 0; op_idx < std::min(sites1_sc.size(), sites2_sc.size()); ++op_idx)
			{
				// get position of the site in the supercell
				const auto [sc1_ok, site1_sc_idx, sc1] = tl2::get_supercell(
					sites1_sc[op_idx], sites_uc, 3, m_eps);
				const auto [sc2_ok, site2_sc_idx, sc2] = tl2::get_supercell(
					sites2_sc[op_idx], sites_uc, 3, m_eps);

				if(!sc1_ok || !sc2_ok)
				{
					std::cerr << "Magdyn error: Could not find supercell for position generated from symop "
						<< op_idx << "." << std::endl;
				}

				ExchangeTerm newterm = term;
				newterm.site1_calc = site1_sc_idx;
				newterm.site2_calc = site2_sc_idx;
				newterm.site1 = GetMagneticSite(newterm.site1_calc).name;
				newterm.site2 = GetMagneticSite(newterm.site2_calc).name;
				newterm.dist_calc = to_3vec<t_vec_real>(sc2 - sc1);
				newterm.dist[0] = tl2::var_to_str(newterm.dist_calc[0], m_prec);
				newterm.dist[1] = tl2::var_to_str(newterm.dist_calc[1], m_prec);
				newterm.dist[2] = tl2::var_to_str(newterm.dist_calc[2], m_prec);

				for(std::uint8_t idx1 = 0; idx1 < 3; ++idx1)
				{
					newterm.dmi[idx1] = tl2::var_to_str(newdmis[op_idx][idx1], m_prec);
					for(std::uint8_t idx2 = 0; idx2 < 3; ++idx2)
						newterm.Jgen[idx1][idx2] = tl2::var_to_str(newJgens[op_idx](idx1, idx2), m_prec);
				}
				newterm.name += "_" + tl2::var_to_str(op_idx + 1, m_prec);

				newterms.emplace_back(std::move(newterm));
			}
		}

		m_exchange_terms = std::move(newterms);
		RemoveDuplicateExchangeTerms();
		CalcSymmetryIndices(symops);
		CalcExchangeTerms();
	}



	/**
	 * generate possible couplings up to a certain distance
	 */
	void GeneratePossibleExchangeTerms(
		t_real dist_max, t_size _sc_max,
		std::optional<t_size> couplings_max)
	{
		if(GetMagneticSitesCount() == 0)
			return;

		// helper struct for finding possible couplings
		struct PossibleCoupling
		{
			// corresponding unit cell position
			t_vec_real pos1_uc{};
			t_vec_real pos2_uc{};

			// magnetic site position in supercell
			t_vec_real sc_vec{};
			t_vec_real pos2_sc{};

			// coordinates in orthogonal lab units
			t_vec_real pos1_uc_lab{};
			t_vec_real pos2_sc_lab{};

			// corresponding unit cell index
			t_size idx1_uc{};
			t_size idx2_uc{};

			// distance between the two magnetic sites
			t_real dist{};
		};

		CalcExternalField();
		CalcMagneticSites();
		CalcExchangeTerms();

		ExchangeTerms newterms;
		std::vector<PossibleCoupling> couplings;

		// generate a list of supercell vectors
		t_real sc_max = t_real(_sc_max);
		for(t_real sc_h = -sc_max; sc_h <= sc_max; sc_h += 1.)
		for(t_real sc_k = -sc_max; sc_k <= sc_max; sc_k += 1.)
		for(t_real sc_l = -sc_max; sc_l <= sc_max; sc_l += 1.)
		{
			// super cell vector
			const t_vec_real sc_vec = tl2::create<t_vec_real>({ sc_h, sc_k, sc_l });

			for(t_size idx1 = 0; idx1 < GetMagneticSitesCount() - 1; ++idx1)
			for(t_size idx2 = idx1 + 1; idx2 < GetMagneticSitesCount(); ++idx2)
			{
				PossibleCoupling coupling;

				coupling.idx1_uc = idx1;
				coupling.idx2_uc = idx2;

				coupling.pos1_uc = GetMagneticSite(idx1).pos_calc;
				coupling.pos2_uc = GetMagneticSite(idx2).pos_calc;

				coupling.sc_vec = sc_vec;
				coupling.pos2_sc = coupling.pos2_uc + sc_vec;

				// transform to lab units for correct distance
				coupling.pos1_uc_lab = m_xtalA * coupling.pos1_uc;
				coupling.pos2_sc_lab = m_xtalA * coupling.pos2_sc;

				coupling.dist = tl2::norm<t_vec_real>(
					coupling.pos2_sc_lab - coupling.pos1_uc_lab);
				if(coupling.dist <= dist_max && coupling.dist > m_eps)
					couplings.emplace_back(std::move(coupling));
			}
		}

		// sort couplings by distance
		std::stable_sort(couplings.begin(), couplings.end(),
			[](const PossibleCoupling& coupling1, const PossibleCoupling& coupling2) -> bool
		{
			return coupling1.dist < coupling2.dist;
		});

		// add couplings to list
		t_size coupling_idx = 0;
		for(const PossibleCoupling& coupling : couplings)
		{
			if(couplings_max && coupling_idx >= *couplings_max)
				break;

			ExchangeTerm newterm{};
			newterm.site1_calc = coupling.idx1_uc;
			newterm.site2_calc = coupling.idx2_uc;
			newterm.site1 = GetMagneticSite(newterm.site1_calc).name;
			newterm.site2 = GetMagneticSite(newterm.site2_calc).name;
			newterm.dist_calc = coupling.sc_vec;
			newterm.dist[0] = tl2::var_to_str(newterm.dist_calc[0], m_prec);
			newterm.dist[1] = tl2::var_to_str(newterm.dist_calc[1], m_prec);
			newterm.dist[2] = tl2::var_to_str(newterm.dist_calc[2], m_prec);
			newterm.length_calc = coupling.dist;
			newterm.J = "0";
			newterm.name += "coupling_" + tl2::var_to_str(coupling_idx + 1, m_prec);
			newterms.emplace_back(std::move(newterm));

			++coupling_idx;
		}

		m_exchange_terms = std::move(newterms);
		RemoveDuplicateExchangeTerms();
		//CalcSymmetryIndices(symops);
		CalcExchangeTerms();
	}



	/**
	 * extend the magnetic structure
	 */
	void ExtendStructure(t_size x_size, t_size y_size, t_size z_size,
		bool remove_duplicates = true, bool flip_spin = true)
	{
		CalcExternalField();
		CalcMagneticSites();

		const t_size num_sites = GetMagneticSitesCount();
		const t_size num_terms = GetExchangeTermsCount();
		m_sites.reserve(num_sites * x_size * y_size * z_size);
		m_exchange_terms.reserve(num_terms * x_size * y_size * z_size);

		// iterate over extended structure
		for(t_size x_idx = 0; x_idx < x_size; ++x_idx)
		for(t_size y_idx = 0; y_idx < y_size; ++y_idx)
		for(t_size z_idx = 0; z_idx < z_size; ++z_idx)
		{
			// ignore sites in original cell
			if(x_idx == 0 && y_idx == 0 && z_idx == 0)
				continue;

			std::string ext_id = "_" + tl2::var_to_str(x_idx + 1, m_prec) +
				"_" + tl2::var_to_str(y_idx + 1, m_prec) +
				"_" + tl2::var_to_str(z_idx + 1, m_prec);

			// duplicate sites
			for(t_size site_idx = 0; site_idx < num_sites; ++site_idx)
			{
				MagneticSite newsite = GetMagneticSite(site_idx);

				newsite.name += ext_id;
				newsite.pos_calc += tl2::create<t_vec_real>({
					t_real(x_idx), t_real(y_idx), t_real(z_idx) });
				newsite.pos[0] = tl2::var_to_str(newsite.pos_calc[0], m_prec);
				newsite.pos[1] = tl2::var_to_str(newsite.pos_calc[1], m_prec);
				newsite.pos[2] = tl2::var_to_str(newsite.pos_calc[2], m_prec);

				if(flip_spin && (x_idx + y_idx + z_idx) % 2 != 0)
				{
					// flip spin
					newsite.spin_dir_calc = -newsite.spin_dir_calc;
					newsite.spin_dir[0] = tl2::var_to_str(newsite.spin_dir_calc[0], m_prec);
					newsite.spin_dir[1] = tl2::var_to_str(newsite.spin_dir_calc[1], m_prec);
					newsite.spin_dir[2] = tl2::var_to_str(newsite.spin_dir_calc[2], m_prec);
				}

				AddMagneticSite(std::move(newsite));
			}

			// duplicate couplings
			for(t_size term_idx = 0; term_idx < num_terms; ++term_idx)
			{
				ExchangeTerm newterm = GetExchangeTerm(term_idx);

				newterm.site1 += ext_id;
				newterm.site2 += ext_id;

				AddExchangeTerm(std::move(newterm));
			}
		}

		if(remove_duplicates)
		{
			RemoveDuplicateMagneticSites();
			RemoveDuplicateExchangeTerms();
		}

		FixExchangeTerms(x_size, y_size, z_size);
		CalcMagneticSites();
		CalcExchangeTerms();
	}



	/**
	 * modify exchange term whose sites point to sc positions that are also available in the uc
	 */
	void FixExchangeTerms(t_size x_size = 1, t_size y_size = 1, t_size z_size = 1)
	{
		for(ExchangeTerm& term : GetExchangeTerms())
		{
			// coupling within uc?
			if(tl2::equals_0<t_vec_real>(term.dist_calc, m_eps))
				continue;

			// find site 2
			const MagneticSite* site2_uc = FindMagneticSite(term.site2);
			if(!site2_uc)
				continue;

			// get site 2 sc vector
			t_vec_real site2_sc = site2_uc->pos_calc + term.dist_calc;

			// fix couplings that are now internal:
			// see if site 2's sc vector is now also available in the uc
			bool fixed_coupling = false;
			for(const MagneticSite& site : GetMagneticSites())
			{
				if(tl2::equals<t_vec_real>(site.pos_calc, site2_sc, m_eps))
				{
					// found the identical site
					term.site2 = site.name;
					term.dist[0] = term.dist[1] = term.dist[2] = "0";
					term.dist_calc = tl2::zero<t_vec_real>(3);
					fixed_coupling = true;
					break;
				}
			}
			if(fixed_coupling)
				continue;

			// fix external couplings
			t_vec_real site2_newsc = site2_sc;

			site2_newsc[0] = std::fmod(site2_newsc[0], t_real(x_size));
			site2_newsc[1] = std::fmod(site2_newsc[1], t_real(y_size));
			site2_newsc[2] = std::fmod(site2_newsc[2], t_real(z_size));

			if(site2_newsc[0] < t_real(0))
				site2_newsc[0] += t_real(x_size);
			if(site2_newsc[1] < t_real(0))
				site2_newsc[1] += t_real(y_size);
			if(site2_newsc[2] < t_real(0))
				site2_newsc[2] += t_real(z_size);

			for(const MagneticSite& site : GetMagneticSites())
			{
				if(tl2::equals<t_vec_real>(site.pos_calc, site2_newsc, m_eps))
				{
					// found the identical site
					term.site2 = site.name;
					term.dist_calc = site2_sc - site2_newsc;
					term.dist[0] = tl2::var_to_str(term.dist_calc[0], m_prec);
					term.dist[1] = tl2::var_to_str(term.dist_calc[1], m_prec);
					term.dist[2] = tl2::var_to_str(term.dist_calc[2], m_prec);
					break;
				}
			}
		}
	}



	/**
	 * remove literal duplicate sites (not symmetry-equivalent ones)
	 */
	void RemoveDuplicateMagneticSites()
	{
		for(auto iter1 = m_sites.begin(); iter1 != m_sites.end(); std::advance(iter1, 1))
		for(auto iter2 = std::next(iter1, 1); iter2 != m_sites.end();)
		{
			if(tl2::equals<t_vec_real>(iter1->pos_calc, iter2->pos_calc, m_eps))
				iter2 = m_sites.erase(iter2);
			else
				std::advance(iter2, 1);
		}
	}



	/**
	 * remove literal duplicate couplings (not symmetry-equivalent ones)
	 */
	void RemoveDuplicateExchangeTerms()
	{
		for(auto iter1 = m_exchange_terms.begin(); iter1 != m_exchange_terms.end(); std::advance(iter1, 1))
		for(auto iter2 = std::next(iter1, 1); iter2 != m_exchange_terms.end();)
		{
			// identical coupling
			bool same_uc = (iter1->site1 == iter2->site1 && iter1->site2 == iter2->site2);
			bool same_sc = tl2::equals<t_vec_real>(iter1->dist_calc, iter2->dist_calc, m_eps);

			// flipped coupling
			bool inv_uc = (iter1->site1 == iter2->site2 && iter1->site2 == iter2->site1);
			bool inv_sc = tl2::equals<t_vec_real>(iter1->dist_calc, -iter2->dist_calc, m_eps);

			if((same_uc && same_sc) || (inv_uc && inv_sc))
				iter2 = m_exchange_terms.erase(iter2);
			else
				std::advance(iter2, 1);
		}
	}



	/**
	 * are two sites equivalent with respect to the given symmetry operators?
	 */
	bool IsSymmetryEquivalent(const MagneticSite& site1, const MagneticSite& site2,
		const std::vector<t_mat_real>& symops) const
	{
		// get symmetry-equivalent positions
		const auto positions = tl2::apply_ops_hom<t_vec_real, t_mat_real, t_real>(
			site1.pos_calc, symops, m_eps);

		for(const auto& pos : positions)
		{
			// symmetry-equivalent site found?
			if(tl2::equals<t_vec_real>(site2.pos_calc, pos, m_eps))
				return true;
		}

		return false;
	}



	/**
	 * are two couplings equivalent with respect to the given symmetry operators?
	 */
	bool IsSymmetryEquivalent(const ExchangeTerm& term1, const ExchangeTerm& term2,
		const std::vector<t_mat_real>& symops) const
	{
		// check if site indices are within bounds
		const t_size N = GetMagneticSitesCount();
		if(term1.site1_calc >= N || term1.site2_calc >= N ||
			term2.site1_calc >= N || term2.site2_calc >= N)
			return false;

		// create unit cell site vectors
		std::vector<t_vec_real> sites_uc = GetMagneticSitePositions(true);

		// super cell distance vector
		t_vec_real dist_sc = to_4vec<t_vec_real>(term1.dist_calc, 0.);

		// generate new (possibly supercell) sites with symop
		auto sites1_sc = tl2::apply_ops_hom<t_vec_real, t_mat_real, t_real>(
			sites_uc[term1.site1_calc], symops, m_eps,
			false /*keep in uc*/, true /*ignore occupied*/,
			true /*return homogeneous*/);
		auto sites2_sc = tl2::apply_ops_hom<t_vec_real, t_mat_real, t_real>(
			sites_uc[term1.site2_calc] + dist_sc, symops, m_eps,
			false /*keep in uc*/, true /*ignore occupied*/,
			true /*return homogeneous*/);

		for(t_size idx = 0; idx < std::min(sites1_sc.size(), sites2_sc.size()); ++idx)
		{
			// get position of the site in the supercell
			const auto [sc1_ok, site1_sc_idx, sc1] = tl2::get_supercell(sites1_sc[idx], sites_uc, 3, m_eps);
			const auto [sc2_ok, site2_sc_idx, sc2] = tl2::get_supercell(sites2_sc[idx], sites_uc, 3, m_eps);
			if(!sc1_ok || !sc2_ok)
				continue;

			// symmetry-equivalent coupling found?
			if(tl2::equals<t_vec_real>(to_3vec<t_vec_real>(sc2 - sc1), term2.dist_calc, m_eps)
				&& site1_sc_idx == term2.site1_calc && site2_sc_idx == term2.site2_calc)
				return true;
		}

		return false;
	}



	/**
	 * assign symmetry group indices to sites and couplings
	 */
	void CalcSymmetryIndices(const std::vector<t_mat_real>& symops)
	{
		// iterate sites
		t_size site_sym_idx_ctr = 0;
		std::vector<const MagneticSite*> seen_sites;

		for(MagneticSite& site : GetMagneticSites())
		{
			const MagneticSite* equivalent_site = nullptr;
			for(const MagneticSite* seen_site : seen_sites)
			{
				// found symmetry-equivalent site
				if(IsSymmetryEquivalent(site, *seen_site, symops))
				{
					equivalent_site = seen_site;
					break;
				}
			}

			if(equivalent_site)
			{
				// assign symmetry index from equivalent site
				site.sym_idx = equivalent_site->sym_idx;
			}
			else
			{
				// assign new symmetry index
				site.sym_idx = ++site_sym_idx_ctr;
				seen_sites.push_back(&site);
			}
		}


		// iterate couplings
		std::vector<const ExchangeTerm*> seen_terms;
		t_size term_sym_idx_ctr = 0;

		for(ExchangeTerm& term : GetExchangeTerms())
		{
			const ExchangeTerm* equivalent_term = nullptr;
			for(const ExchangeTerm* seen_term : seen_terms)
			{
				// found symmetry-equivalent coupling
				if(IsSymmetryEquivalent(term, *seen_term, symops))
				{
					equivalent_term = seen_term;
					break;
				}
			}

			if(equivalent_term)
			{
				// assign symmetry index from equivalent coupling
				term.sym_idx = equivalent_term->sym_idx;
			}
			else
			{
				// assign new symmetry index
				term.sym_idx = ++term_sym_idx_ctr;
				seen_terms.push_back(&term);
			}
		}
	}
	// --------------------------------------------------------------------



	// --------------------------------------------------------------------
	// calculation functions
	// --------------------------------------------------------------------
	/**
	 * calculate the rotation matrix for the external field
	 */
	void CalcExternalField()
	{
		bool use_field =
			(!tl2::equals_0<t_real>(m_field.mag, m_eps) || m_field.align_spins)
			&& m_field.dir.size() == 3;

		if(use_field)
		{
			// rotate field to [001] direction
			m_rot_field = tl2::convert<t_mat>(
				tl2::trans<t_mat_real>(
					tl2::rotation<t_mat_real, t_vec_real>(
						-m_field.dir, m_zdir, &m_rotaxis, m_eps)));

#ifdef __TLIBS2_MAGDYN_DEBUG_OUTPUT__
			std::cout << "Field rotation from:\n";
			tl2::niceprint(std::cout, -m_field.dir, 1e-4, 4);
			std::cout << "\nto:\n";
			tl2::niceprint(std::cout, m_zdir, 1e-4, 4);
			std::cout << "\nmatrix:\n";
			tl2::niceprint(std::cout, m_rot_field, 1e-4, 4);
			std::cout << std::endl;
#endif
		}
	}



	/**
	 * calculate the spin rotation trafo for the magnetic sites
	 * and parse any given expressions
	 */
	void CalcMagneticSite(MagneticSite& site)
	{
		try
		{
			tl2::ExprParser<t_cplx> parser = GetExprParser();

			// spin direction and orthogonal plane
			bool has_explicit_trafo = true;

			// defaults
			site.pos_calc = tl2::zero<t_vec_real>(3);
			site.spin_dir_calc = tl2::zero<t_vec_real>(3);
			site.trafo_z_calc = tl2::zero<t_vec>(3);
			site.trafo_plane_calc = tl2::zero<t_vec>(3);
			site.trafo_plane_conj_calc = tl2::zero<t_vec>(3);
			if(site.g_e.size1() == 0 || site.g_e.size2() == 0)
				site.g_e = tl2::g_e<t_real> * tl2::unit<t_mat>(3);

			// spin magnitude
			if(parser.parse_noexcept(site.spin_mag))
			{
				site.spin_mag_calc = parser.eval_noexcept().real();
			}
			else
			{
				std::cerr << "Magdyn error: Parsing spin magnitude \""
					<< site.spin_mag << "\""
					<< " for site \"" << site.name << "\""
					<< "." << std::endl;
			}


			for(std::uint8_t idx = 0; idx < 3; ++idx)
			{
				// position
				if(site.pos[idx] != "")
				{
					if(parser.parse_noexcept(site.pos[idx]))
					{
						site.pos_calc[idx] = parser.eval_noexcept().real();
					}
					else
					{
						std::cerr << "Magdyn error: Parsing position \""
							<< site.pos[idx] << "\""
							<< " for site \"" << site.name << "\""
							<< " and component " << idx
							<< "." << std::endl;
					}
				}

				// spin direction
				if(site.spin_dir[idx] != "")
				{
					if(parser.parse_noexcept(site.spin_dir[idx]))
					{
						site.spin_dir_calc[idx] = parser.eval_noexcept().real();
					}
					else
					{
						std::cerr << "Magdyn error: Parsing spin direction \""
							<< site.spin_dir[idx] << "\""
							<< " for site \"" << site.name << "\""
							<< " and component " << idx
							<< "." << std::endl;
					}
				}

				// orthogonal spin direction
				if(site.spin_ortho[idx] != "")
				{
					if(parser.parse_noexcept(site.spin_ortho[idx]))
					{
						site.trafo_plane_calc[idx] = parser.eval_noexcept();
						site.trafo_plane_conj_calc[idx] = std::conj(site.trafo_plane_calc[idx]);
					}
					else
					{
						has_explicit_trafo = false;

						std::cerr << "Magdyn error: Parsing spin orthogonal plane \""
							<< site.spin_ortho[idx] << "\""
							<< " for site \"" << site.name << "\""
							<< " and component " << idx
							<< "." << std::endl;
					}
				}
				else
				{
					has_explicit_trafo = false;
				}
			}

			// spin rotation of equation (9) from (Toth 2015)
			if(m_field.align_spins)
			{
				std::tie(site.trafo_plane_calc, site.trafo_z_calc) =
					rot_to_trafo(m_rot_field);
			}
			else
			{
				if(!has_explicit_trafo)
				{
					// calculate u and v from the spin rotation
					std::tie(site.trafo_plane_calc, site.trafo_z_calc) =
						spin_to_trafo(site.spin_dir_calc);
				}

				// TODO: normalise the v vector as well as the real and imaginary u vectors
				// in case they are explicitly given

#ifdef __TLIBS2_MAGDYN_DEBUG_OUTPUT__
				std::cout << "Site " << site.name << " u = "
					<< site.trafo_plane_calc[0] << " " << site.trafo_plane_calc[1] << " " << site.trafo_plane_calc[2]
					<< std::endl;
				std::cout << "Site " << site.name << " v = "
					<< site.trafo_z_calc[0] << " " << site.trafo_z_calc[1] << " " << site.trafo_z_calc[2]
					<< std::endl;
#endif
			}

			site.trafo_plane_conj_calc = tl2::conj(site.trafo_plane_calc);

			// multiply g factor
			site.ge_trafo_z_calc = site.g_e * site.trafo_z_calc;
			site.ge_trafo_plane_calc = site.g_e * site.trafo_plane_calc;
			site.ge_trafo_plane_conj_calc = site.g_e * site.trafo_plane_conj_calc;

		}
		catch(const std::exception& ex)
		{
			std::cerr << "Magdyn error: Calculating site \"" << site.name << "\"."
				<< " Reason: " << ex.what()
				<< std::endl;
		}
	}



	/**
	 * calculate the spin rotation trafo for the magnetic sites
	 * and parse any given expressions
	 */
	void CalcMagneticSites()
	{
		for(MagneticSite& site : GetMagneticSites())
			CalcMagneticSite(site);
	}



	/**
	 * parse the exchange term expressions and calculate all properties
	 */
	void CalcExchangeTerm(ExchangeTerm& term)
	{
		try
		{
			tl2::ExprParser<t_cplx> parser = GetExprParser();

			// defaults
			term.dist_calc = tl2::zero<t_vec_real>(3);  // distance
			term.dmi_calc = tl2::zero<t_vec>(3);        // dmi interaction
			term.Jgen_calc = tl2::zero<t_mat>(3, 3);    // general exchange interaction

			// get site indices
			term.site1_calc = GetMagneticSiteIndex(term.site1);
			term.site2_calc = GetMagneticSiteIndex(term.site2);

#ifdef __TLIBS2_MAGDYN_DEBUG_OUTPUT__
			std::cout << "Coupling: "
				<< term.name << ": " << term.site1 << " (" << term.site1_calc << ") -> "
				<< term.site2 << " (" << term.site2_calc << ")." << std::endl;
#endif

			if(term.site1_calc >= GetMagneticSitesCount())
			{
				std::cerr << "Magdyn error: Unknown site 1 name \"" << term.site1 << "\"."
					<< " in coupling \"" << term.name << "\"."
					<< std::endl;
				return;
			}
			if(term.site2_calc >= GetMagneticSitesCount())
			{
				std::cerr << "Magdyn error: Unknown site 2 name \"" << term.site2 << "\"."
					<< " in coupling \"" << term.name << "\"."
					<< std::endl;
				return;
			}

			// symmetric interaction
			if(term.J == "")
			{
				term.J_calc = t_real(0);
			}
			else if(parser.parse_noexcept(term.J))
			{
				term.J_calc = parser.eval_noexcept();
			}
			else
			{
				std::cerr << "Magdyn error: Parsing J term \""
					<< term.J << "\"." << std::endl;
			}

			for(std::uint8_t i = 0; i < 3; ++i)
			{
				// distance
				if(term.dist[i] != "")
				{
					if(parser.parse_noexcept(term.dist[i]))
					{
						term.dist_calc[i] = parser.eval_noexcept().real();
					}
					else
					{
						std::cerr << "Magdyn error: Parsing distance term \""
							<< term.dist[i]
							<< "\" (index " << i << ")"
							<< "." << std::endl;
					}
				}

				// dmi
				if(term.dmi[i] != "")
				{
					if(parser.parse_noexcept(term.dmi[i]))
					{
						term.dmi_calc[i] = parser.eval_noexcept();
					}
					else
					{
						std::cerr << "Magdyn error: Parsing DMI term \""
							<< term.dmi[i]
							<< "\" (index " << i << ")"
							<< "." << std::endl;
					}
				}

				// general exchange interaction
				for(std::uint8_t j = 0; j < 3; ++j)
				{
					if(term.Jgen[i][j] == "")
						continue;

					if(parser.parse_noexcept(term.Jgen[i][j]))
					{
						term.Jgen_calc(i, j) = parser.eval_noexcept();
					}
					else
					{
						std::cerr << "Magdyn error: Parsing general term \""
							<< term.Jgen[i][j]
							<< "\" (indices " << i << ", " << j << ")"
							<< "." << std::endl;
					}
				}
			}

			const t_vec_real& pos1_uc = GetMagneticSite(term.site1_calc).pos_calc;
			const t_vec_real& pos2_uc = GetMagneticSite(term.site2_calc).pos_calc;
			t_vec_real pos2_sc = pos2_uc + term.dist_calc;

			// transform to lab units for correct distance
			t_vec_real pos1_uc_lab = m_xtalA * pos1_uc;
			t_vec_real pos2_sc_lab = m_xtalA * pos2_sc;

			term.length_calc = tl2::norm<t_vec_real>(pos2_sc_lab - pos1_uc_lab);
		}
		catch(const std::exception& ex)
		{
			std::cerr << "Magdyn error: Calculating coupling \"" << term.name << "\"."
				<< " Reason: " << ex.what() << "." << std::endl;
		}
	}



	/**
	 * parse all exchange term expressions and calculate all properties
	 */
	void CalcExchangeTerms()
	{
		for(ExchangeTerm& term : GetExchangeTerms())
			CalcExchangeTerm(term);
	}



	/**
	 * calculate the real-space interaction matrix J of
	 * equations (10) - (13) from (Toth 2015)
	 */
	t_mat CalcRealJ(const ExchangeTerm& term) const
	{
		// symmetric part of the exchange interaction matrix, see (Toth 2015) p. 2
		t_mat J = tl2::diag<t_mat>(
			tl2::create<t_vec>({ term.J_calc, term.J_calc, term.J_calc }));

		// dmi as anti-symmetric part of interaction matrix
		// using a cross product matrix, see (Toth 2015) p. 2
		if(term.dmi_calc.size() == 3)
			J += tl2::skewsymmetric<t_mat, t_vec>(-term.dmi_calc);

		// general J matrix
		if(term.Jgen_calc.size1() == 3 && term.Jgen_calc.size2() == 3)
			J += term.Jgen_calc;

		// incommensurate case: rotation wrt magnetic unit cell
		// equations (21), (6), (2) as well as section 10 from (Toth 2015)
		if(IsIncommensurate())
		{
			const t_real rot_UC_angle =
				s_twopi * tl2::inner<t_vec_real>(m_ordering, term.dist_calc);

			if(!tl2::equals_0<t_real>(rot_UC_angle, m_eps))
			{
				t_mat rot_UC = tl2::convert<t_mat>(
					tl2::rotation<t_mat_real, t_vec_real>(
						m_rotaxis, rot_UC_angle));
				J = J * rot_UC;

#ifdef __TLIBS2_MAGDYN_DEBUG_OUTPUT__
				std::cout << "Coupling rot_UC = " << term.name << ":\n";
				tl2::niceprint(std::cout, rot_UC, 1e-4, 4);
#endif
			}
		}

		return J;
	}



	/**
	 * calculate the reciprocal interaction matrices J(Q) and J(-Q) of
	 * equations (12) and (14) from (Toth 2015)
	 */
	std::tuple<t_Jmap, t_Jmap> CalcReciprocalJs(const t_vec_real& Qvec) const
	{
		t_Jmap J_Q, J_Q0;

		// no (or not valid) exchange terms given
		if(GetExchangeTermsCount() == 0)
			return std::make_tuple(J_Q, J_Q0);

		// iterate couplings to pre-calculate corresponding J matrices
		for(const ExchangeTerm& term : GetExchangeTerms())
		{
			// insert or add an exchange matrix at the given site indices
			auto insert_or_add = [](t_Jmap& J, const t_indices& indices, const t_mat& J33)
			{
				if(auto iter = J.find(indices); iter != J.end())
					iter->second += J33;
				else
					J.emplace(std::move(std::make_pair(indices, J33)));
			};

			if(!CheckMagneticSite(term.site1_calc) || !CheckMagneticSite(term.site2_calc))
				continue;

			const t_indices indices = std::make_pair(term.site1_calc, term.site2_calc);
			const t_indices indices_t = std::make_pair(term.site2_calc, term.site1_calc);

			const t_mat J = CalcRealJ(term);
			if(J.size1() == 0 || J.size2() == 0)
				continue;
			const t_mat J_T = tl2::trans(J);

			// get J in reciprocal space by fourier trafo
			// equations (14), (12), (11), and (52) from (Toth 2015)
			const t_cplx phase = m_phase_sign * s_imag * s_twopi *
				tl2::inner<t_vec_real>(term.dist_calc, Qvec);

			insert_or_add(J_Q, indices, J * std::exp(phase));
			insert_or_add(J_Q, indices_t, J_T * std::exp(-phase));

			insert_or_add(J_Q0, indices, J);
			insert_or_add(J_Q0, indices_t, J_T);
		}  // end of iteration over couplings

		return std::make_tuple(J_Q, J_Q0);
	}



	/**
	 * get the hamiltonian at the given momentum
	 * @note implements the formalism given by (Toth 2015)
	 * @note a first version for a simplified ferromagnetic dispersion was based on (Heinsdorf 2021)
	 */
	t_mat CalcHamiltonian(const t_vec_real& Qvec) const
	{
		const t_size N = GetMagneticSitesCount();
		if(N == 0)
			return t_mat{};

		// build the interaction matrices J(Q) and J(-Q) of
		// equations (12) and (14) from (Toth 2015)
		const auto [J_Q, J_Q0] = CalcReciprocalJs(Qvec);

		// create the hamiltonian of equation (25) and (26) from (Toth 2015)
		t_mat H00     = tl2::zero<t_mat>(N, N);
		t_mat H00c_mQ = tl2::zero<t_mat>(N, N);  // H00*(-Q)
		t_mat H0N     = tl2::zero<t_mat>(N, N);

		bool use_field = !tl2::equals_0<t_real>(m_field.mag, m_eps)
			&& m_field.dir.size() == 3;

		// iterate magnetic sites
		for(t_size i = 0; i < N; ++i)
		{
			const MagneticSite& s_i = GetMagneticSite(i);

			// get the pre-calculated u and v vectors for the commensurate case
			const t_vec& u_i  = s_i.trafo_plane_calc;
			const t_vec& uc_i = s_i.trafo_plane_conj_calc;  // u*_i
			const t_vec& v_i  = s_i.trafo_z_calc;

			for(t_size j = 0; j < N; ++j)
			{
				const MagneticSite& s_j = GetMagneticSite(j);

				// get the pre-calculated u and v vectors for the commensurate case
				const t_vec& u_j  = s_j.trafo_plane_calc;
				const t_vec& uc_j = s_j.trafo_plane_conj_calc;  // u*_j
				const t_vec& v_j  = s_j.trafo_z_calc;

				// get the pre-calculated exchange matrices for the (i, j) coupling
				const t_indices indices_ij = std::make_pair(i, j);
				const t_mat* J_Q33 = nullptr;
				const t_mat* J_Q033 = nullptr;
				if(auto iter = J_Q.find(indices_ij); iter != J_Q.end())
					J_Q33 = &iter->second;
				if(auto iter = J_Q0.find(indices_ij); iter != J_Q0.end())
					J_Q033 = &iter->second;

				if(J_Q33)
				{
					const t_real S_mag = 0.5 * std::sqrt(s_i.spin_mag_calc * s_j.spin_mag_calc);

					// equation (26) from (Toth 2015)
					H00(i, j)     += S_mag * tl2::inner_noconj<t_vec>(u_i,  (*J_Q33) * uc_j);
					H00c_mQ(i, j) += S_mag * tl2::inner_noconj<t_vec>(uc_i, (*J_Q33) * u_j);
					H0N(i, j)     += S_mag * tl2::inner_noconj<t_vec>(u_i,  (*J_Q33) * u_j);
				}

				if(J_Q033)
				{
					// equation (26) from (Toth 2015)
					t_cplx c = s_j.spin_mag_calc * tl2::inner_noconj<t_vec>(v_i, (*J_Q033) * v_j);
					H00(i, i)     -= c;
					H00c_mQ(i, i) -= c;
				}
			}  // end of iteration over j sites

			// include external field, equation (28) from (Toth 2015)
			if(use_field)
			{
				const t_vec field = tl2::convert<t_vec>(-m_field.dir) * m_field.mag;
				const t_vec gv    = s_i.g_e * v_i;
				const t_cplx Bgv  = tl2::inner_noconj<t_vec>(field, gv);

				// bohr magneton in [meV/T]
				constexpr const t_real muB = tl2::mu_B<t_real>
					/ tl2::meV<t_real> * tl2::tesla<t_real>;

				H00(i, i)     -= muB * Bgv;
				H00c_mQ(i, i) -= std::conj(muB * Bgv);
			}
		}  // end of iteration over i sites

		// equation (25) from (Toth 2015)
		t_mat H = tl2::create<t_mat>(N*2, N*2);
		tl2::set_submat(H, H00,            0, 0);
		tl2::set_submat(H, H0N,            0, N);
		tl2::set_submat(H, tl2::herm(H0N), N, 0);
		tl2::set_submat(H, H00c_mQ,        N, N);

		return H;
	}



	/**
	 * get the energies from a hamiltonian
	 * @note implements the formalism given by (Toth 2015)
	 */
	EnergiesAndWeights CalcEnergiesFromHamiltonian(
		t_mat _H, const t_vec_real& Qvec,
		bool only_energies = false) const
	{
		using namespace tl2_ops;
		const t_size N = GetMagneticSitesCount();
		if(N == 0 || _H.size1() == 0 || _H.size2() == 0)
			return {};

		// equation (30) from (Toth 2015)
		t_mat g_sign = tl2::unit<t_mat>(N*2, N*2);
		for(t_size i = N; i < 2*N; ++i)
			g_sign(i, i) = -1.;

		// equation (31) from (Toth 2015)
		t_mat C_mat;
		t_size chol_try = 0;
		for(; chol_try < m_tries_chol; ++chol_try)
		{
			const auto [chol_ok, _C] = tl2_la::chol<t_mat>(_H);

			if(chol_ok)
			{
				C_mat = std::move(_C);
				break;
			}
			else
			{
				if(chol_try >= m_tries_chol - 1)
				{
					std::cerr << "Magdyn warning: Cholesky decomposition failed at Q = "
						<< Qvec << "." << std::endl;
					C_mat = std::move(_C);
					break;
				}

				// try forcing the hamilton to be positive definite
				for(t_size i = 0; i < 2*N; ++i)
					_H(i, i) += m_delta_chol;
			}
		}

		if(m_perform_checks && chol_try > 0)
		{
			std::cerr << "Magdyn warning: Needed " << chol_try
				<< " correction(s) for Cholesky decomposition at Q = "
				<< Qvec << "." << std::endl;
		}

		if(C_mat.size1() == 0 || C_mat.size2() == 0)
		{
			std::cerr << "Magdyn error: Invalid Cholesky decomposition at Q = "
				<< Qvec << "." << std::endl;
			return {};
		}

		// see p. 5 in (Toth 2015)
		const t_mat H_mat = C_mat * g_sign * tl2::herm<t_mat>(C_mat);

		const bool is_herm = tl2::is_symm_or_herm<t_mat, t_real>(H_mat, m_eps);
		if(m_perform_checks && !is_herm)
		{
			std::cerr << "Magdyn warning: Hamiltonian is not hermitian at Q = "
				<< Qvec << "." << std::endl;
		}

		// eigenvalues of the hamiltonian correspond to the energies
		// eigenvectors correspond to the spectral weights
		const auto [evecs_ok, evals, evecs] =
			tl2_la::eigenvec<t_mat, t_vec, t_cplx, t_real>(
				H_mat, only_energies, is_herm, true);
		if(!evecs_ok)
		{
			std::cerr << "Magdyn warning: Eigensystem calculation failed at Q = "
				<< Qvec << "." << std::endl;
		}


		EnergiesAndWeights energies_and_correlations{};
		energies_and_correlations.reserve(evals.size());

		// register energies
		for(const auto& eval : evals)
		{
			const EnergyAndWeight EandS { .E = eval.real(), };
			energies_and_correlations.emplace_back(std::move(EandS));
		}

		// weight factors
		if(!only_energies)
		{
			CalcCorrelationsFromHamiltonian(energies_and_correlations,
				H_mat, C_mat, g_sign, Qvec, evecs);
		}

		return energies_and_correlations;
	}



	/**
	 * get the dynamical structure factor from a hamiltonian
	 * @note implements the formalism given by (Toth 2015)
	 */
	void CalcCorrelationsFromHamiltonian(EnergiesAndWeights& energies_and_correlations,
		const t_mat& H_mat, const t_mat& C_mat, const t_mat& g_sign,
		const t_vec_real& Qvec, const std::vector<t_vec>& evecs) const
	{
		const t_size N = GetMagneticSitesCount();
		if(N == 0)
			return;

		// get the sorting of the energies
		const std::vector<t_size> sorting = tl2::get_perm(
			energies_and_correlations.size(),
			[&energies_and_correlations](t_size idx1, t_size idx2) -> bool
		{
			return energies_and_correlations[idx1].E >=
				energies_and_correlations[idx2].E;
		});

		const t_mat evec_mat = tl2::create<t_mat>(tl2::reorder(evecs, sorting));
		const t_mat evec_mat_herm = tl2::herm(evec_mat);

		// equation (32) from (Toth 2015)
		const t_mat L_mat = evec_mat_herm * H_mat * evec_mat;  // energies
		t_mat E_sqrt = g_sign * L_mat;                         // abs. energies
		for(t_size i = 0; i < E_sqrt.size1(); ++i)
			E_sqrt(i, i) = std::sqrt(E_sqrt(i, i));            // sqrt. of abs. energies

		// re-create energies, to be consistent with the weights
		energies_and_correlations.clear();
		for(t_size i = 0; i < L_mat.size1(); ++i)
		{
			const EnergyAndWeight EandS
			{
				.E = L_mat(i, i).real(),
				.S = tl2::zero<t_mat>(3, 3),
				.S_perp = tl2::zero<t_mat>(3, 3),
			};

			energies_and_correlations.emplace_back(std::move(EandS));
		}

		const auto [C_inv, inv_ok] = tl2::inv(C_mat);
		if(!inv_ok)
		{
			using namespace tl2_ops;
			std::cerr << "Magdyn warning: Inversion failed at Q = "
				<< Qvec << "." << std::endl;
		}

		// equation (34) from (Toth 2015)
		const t_mat trafo = C_inv * evec_mat * E_sqrt;
		const t_mat trafo_herm = tl2::herm(trafo);

#ifdef __TLIBS2_MAGDYN_DEBUG_OUTPUT__
		t_mat D_mat = trafo_herm * H_mat * trafo;
		std::cout << "D = \n";
		tl2::niceprint(std::cout, D_mat, 1e-4, 4);
		std::cout << "\nE = \n";
		tl2::niceprint(std::cout, E_sqrt, 1e-4, 4);
		std::cout << "\nL = \n";
		tl2::niceprint(std::cout, L_mat, 1e-4, 4);
		std::cout << std::endl;
#endif

		// building the spin correlation functions of equation (47) from (Toth 2015)
		for(std::uint8_t x_idx = 0; x_idx < 3; ++x_idx)
		for(std::uint8_t y_idx = 0; y_idx < 3; ++y_idx)
		{
			// equations (44) from (Toth 2015)
			t_mat M0N = tl2::create<t_mat>(N, N);
			t_mat M00 = tl2::create<t_mat>(N, N);
			t_mat MNN = tl2::create<t_mat>(N, N);
			t_mat MN0 = tl2::create<t_mat>(N, N);

			for(t_size i = 0; i < N; ++i)
			for(t_size j = 0; j < N; ++j)
			{
				// get the sites
				const MagneticSite& s_i = GetMagneticSite(i);
				const MagneticSite& s_j = GetMagneticSite(j);

				// get the pre-calculated u vectors
				const t_vec& u_i = s_i.ge_trafo_plane_calc;
				const t_vec& u_j = s_j.ge_trafo_plane_calc;
				const t_vec& uc_i = s_i.ge_trafo_plane_conj_calc;
				const t_vec& uc_j = s_j.ge_trafo_plane_conj_calc;

				// pre-factors of equation (44) from (Toth 2015)
				const t_real S_mag = std::sqrt(s_i.spin_mag_calc * s_j.spin_mag_calc);
				const t_cplx phase = std::exp(-m_phase_sign * s_imag * s_twopi *
					tl2::inner<t_vec_real>(s_j.pos_calc - s_i.pos_calc, Qvec));

				// matrix elements of equation (44) from (Toth 2015)
				M0N(i, j) = phase * S_mag * u_i[x_idx]  * u_j[y_idx];
				M00(i, j) = phase * S_mag * u_i[x_idx]  * uc_j[y_idx];
				MNN(i, j) = phase * S_mag * uc_i[x_idx] * u_j[y_idx];
				MN0(i, j) = phase * S_mag * uc_i[x_idx] * uc_j[y_idx];
			} // end of iteration over sites

			// equation (47) from (Toth 2015)
			t_mat M = tl2::create<t_mat>(N*2, N*2);
			tl2::set_submat(M, M0N, 0, N);
			tl2::set_submat(M, M00, 0, 0);
			tl2::set_submat(M, MNN, N, N);
			tl2::set_submat(M, MN0, N, 0);

			const t_mat M_trafo = trafo_herm * M * trafo;

#ifdef __TLIBS2_MAGDYN_DEBUG_OUTPUT__
			std::cout << "M_trafo for x=" << x_idx << ", y=" << y_idx << ":\n";
			tl2::niceprint(std::cout, M_trafo, 1e-4, 4);
			std::cout << std::endl;
#endif

			for(t_size i = 0; i < energies_and_correlations.size(); ++i)
			{
				(energies_and_correlations[i].S)(x_idx, y_idx) +=
					M_trafo(i, i) / t_real(2*N);
			}
		} // end of coordinate iteration
	}



	/**
	 * applies projectors, form and weight factors to get neutron intensities
	 * @note implements the formalism given by (Toth 2015)
	 */
	void CalcIntensities(const t_vec_real& Q_rlu, EnergiesAndWeights&
		energies_and_correlations) const
	{
		using namespace tl2_ops;
		tl2::ExprParser<t_cplx> magffact = m_magffact;

		for(EnergyAndWeight& E_and_S : energies_and_correlations)
		{
			// apply bose factor
			if(m_temperature >= 0.)
				E_and_S.S *= tl2::bose_cutoff(E_and_S.E, m_temperature, m_bose_cutoff);

			// apply form factor
			if(m_magffact_formula != "")
			{
				// get |Q| in units of A^(-1)
				t_vec_real Q_invA = m_xtalB * Q_rlu;
				t_real Q_abs = tl2::norm<t_vec_real>(Q_invA);

				// evaluate form factor expression
				magffact.register_var("Q", Q_abs);
				t_real ffact = magffact.eval_noexcept().real();
				E_and_S.S *= ffact;
			}

			// apply orthogonal projector for magnetic neutron scattering,
			// see (Shirane 2002), p. 37, equation (2.64)
			t_mat proj_neutron = tl2::ortho_projector<t_mat, t_vec>(Q_rlu, false);
			E_and_S.S_perp = proj_neutron * E_and_S.S * proj_neutron;

			// weights
			E_and_S.S_sum       = tl2::trace<t_mat>(E_and_S.S);
			E_and_S.S_perp_sum  = tl2::trace<t_mat>(E_and_S.S_perp);
			E_and_S.weight_full = std::abs(E_and_S.S_sum.real());
			E_and_S.weight      = std::abs(E_and_S.S_perp_sum.real());

			// TODO: polarisation via blume-maleev equation
			/*static const t_vec_real h_rlu = tl2::create<t_vec_real>({ 1., 0., 0. });
			static const t_vec_real l_rlu = tl2::create<t_vec_real>({ 0., 0., 1. });
			t_vec_real h_lab = m_xtalUB * h_rlu;
			t_vec_real l_lab = m_xtalUB * l_rlu;
			t_vec_real Q_lab = m_xtalUB * Q_rlu;
			t_mat_real rotQ = tl2::rotation<t_mat_real>(h_lab, Q_lab, &l_lab, m_eps, true);
			t_mat_real rotQ_hkl = m_xtalUBinv * rotQ * m_xtalUB;
			const auto [rotQ_hkl_inv, rotQ_hkl_inv_ok] = tl2::inv(rotQ_hkl);
			if(!rotQ_hkl_inv_ok)
				std::cerr << "Magdyn error: Cannot invert Q rotation matrix."
					<< std::endl;

			t_mat rotQ_hkl_cplx = tl2::convert<t_mat, t_mat_real>(rotQ_hkl);
			t_mat rotQ_hkl_inv_cplx = tl2::convert<t_mat, t_mat_real>(rotQ_hkl_inv);

			t_mat S_perp_pol = rotQ_hkl_cplx * E_and_S.S_perp * rotQ_hkl_inv_cplx;
			t_mat S_pol = rotQ_hkl_cplx * E_and_S.S * rotQ_hkl_inv_cplx;*/
		}
	}



	/**
	 * unite degenerate energies and their corresponding eigenstates
	 */
	EnergiesAndWeights UniteEnergies(const EnergiesAndWeights&
		energies_and_correlations) const
	{
		EnergiesAndWeights new_energies_and_correlations{};
		new_energies_and_correlations.reserve(energies_and_correlations.size());

		for(const EnergyAndWeight& curState : energies_and_correlations)
		{
			const t_real curE = curState.E;

			auto iter = std::find_if(
				new_energies_and_correlations.begin(),
				new_energies_and_correlations.end(),
				[curE, this](const EnergyAndWeight& E_and_S) -> bool
			{
				t_real E = E_and_S.E;
				return tl2::equals<t_real>(E, curE, m_eps);
			});

			if(iter == new_energies_and_correlations.end())
			{
				// energy not yet seen
				new_energies_and_correlations.push_back(curState);
			}
			else
			{
				// energy already seen: add correlation matrices and weights
				iter->S           += curState.S;
				iter->S_perp      += curState.S_perp;
				iter->S_sum       += curState.S_sum;
				iter->S_perp_sum  += curState.S_perp_sum;
				iter->weight      += curState.weight;
				iter->weight_full += curState.weight_full;
			}
		}

		return new_energies_and_correlations;
	}



	/**
	 * get the energies and the spin-correlation at the given momentum
	 * (also calculates incommensurate contributions and applies weight factors)
	 * @note implements the formalism given by (Toth 2015)
	 */
	EnergiesAndWeights CalcEnergies(const t_vec_real& Q_rlu,
		bool only_energies = false) const
	{
		EnergiesAndWeights EandWs;

		if(m_calc_H)
		{
			const t_mat H = CalcHamiltonian(Q_rlu);
			EandWs = CalcEnergiesFromHamiltonian(H, Q_rlu, only_energies);
		}

		if(IsIncommensurate())
		{
			// equations (39) and (40) from (Toth 2015)
			const t_mat proj_norm = tl2::convert<t_mat>(
				tl2::projector<t_mat_real, t_vec_real>(m_rotaxis, true));

			t_mat rot_incomm = tl2::unit<t_mat>(3);
			rot_incomm -= s_imag * m_phase_sign * tl2::skewsymmetric<t_mat, t_vec>(m_rotaxis);
			rot_incomm -= proj_norm;
			rot_incomm *= 0.5;

			const t_mat rot_incomm_conj = tl2::conj(rot_incomm);

			EnergiesAndWeights EandWs_p, EandWs_m;

			if(m_calc_Hp)
			{
				const t_mat H_p = CalcHamiltonian(Q_rlu + m_ordering);
				EandWs_p = CalcEnergiesFromHamiltonian(
					H_p, Q_rlu + m_ordering, only_energies);
			}

			if(m_calc_Hm)
			{
				const t_mat H_m = CalcHamiltonian(Q_rlu - m_ordering);
				EandWs_m = CalcEnergiesFromHamiltonian(
					H_m, Q_rlu - m_ordering, only_energies);
			}

			if(!only_energies)
			{
				// formula 40 from (Toth 2015)
				for(EnergyAndWeight& EandW : EandWs)
					EandW.S = EandW.S * proj_norm;
				for(EnergyAndWeight& EandW : EandWs_p)
					EandW.S = EandW.S * rot_incomm;
				for(EnergyAndWeight& EandW : EandWs_m)
					EandW.S = EandW.S * rot_incomm_conj;
			}

			// unite energies and weights
			EandWs.reserve(EandWs.size() + EandWs_p.size() + EandWs_m.size());
			for(EnergyAndWeight& EandW : EandWs_p)
				EandWs.emplace_back(std::move(EandW));
			for(EnergyAndWeight& EandW : EandWs_m)
				EandWs.emplace_back(std::move(EandW));
		}

		if(!only_energies)
			CalcIntensities(Q_rlu, EandWs);

		if(m_unite_degenerate_energies)
			EandWs = UniteEnergies(EandWs);

		if(!only_energies)
			CheckImagWeights(Q_rlu, EandWs);

		return EandWs;
	}



	EnergiesAndWeights CalcEnergies(t_real h, t_real k, t_real l,
		bool only_energies = false) const
	{
		// momentum transfer
		const t_vec_real Qvec = tl2::create<t_vec_real>({ h, k, l });
		return CalcEnergies(Qvec, only_energies);
	}



	/**
	 * generates the dispersion plot along the given q path
	 */
	SofQEs CalcDispersion(t_real h_start, t_real k_start, t_real l_start,
		t_real h_end, t_real k_end, t_real l_end,
		t_size num_qs = 128, t_size num_threads = 4) const
	{
		// thread pool and tasks
		using t_pool = boost::asio::thread_pool;
		using t_task = std::packaged_task<SofQE()>;
		using t_taskptr = std::shared_ptr<t_task>;

		t_pool pool{num_threads};
		std::vector<t_taskptr> tasks;
		tasks.reserve(num_qs);

		// calculate dispersion
		for(t_size i = 0; i < num_qs; ++i)
		{
			auto task = [this, i, num_qs,
				h_start, k_start, l_start,
				h_end, k_end, l_end]() -> SofQE
			{
				// get Q
				const t_real h = std::lerp(h_start, h_end, t_real(i) / t_real(num_qs - 1));
				const t_real k = std::lerp(k_start, k_end, t_real(i) / t_real(num_qs - 1));
				const t_real l = std::lerp(l_start, l_end, t_real(i) / t_real(num_qs - 1));

				// get E and S(Q, E) for this Q
				EnergiesAndWeights E_and_S = CalcEnergies(h, k, l, false);

				SofQE result
				{
					.h = h, .k = k, .l = l,
					.E_and_S = E_and_S
				};

				return result;
			};

			t_taskptr taskptr = std::make_shared<t_task>(task);
			tasks.push_back(taskptr);
			boost::asio::post(pool, [taskptr]() { (*taskptr)(); });
		}

		// collect results
		SofQEs results;
		results.reserve(tasks.size());

		for(auto& task : tasks)
		{
			const SofQE& result = task->get_future().get();
			results.push_back(result);
		}

		return results;
	}



	/**
	 * get the energy minimum
	 * @note a first version for a simplified ferromagnetic dispersion was based on (Heinsdorf 2021)
	 */
	t_real CalcMinimumEnergy() const
	{
		// energies at (000)
		const auto energies_and_correlations = CalcEnergies(0., 0., 0., true);

		// get minimum
		const auto min_iter = std::min_element(
			energies_and_correlations.begin(), energies_and_correlations.end(),
			[](const EnergyAndWeight& E_and_S_1, const EnergyAndWeight& E_and_S_2) -> bool
		{
			return std::abs(E_and_S_1.E) < std::abs(E_and_S_2.E);
		});

		if(min_iter == energies_and_correlations.end())
			return 0.;
		return min_iter->E;
	}



	/**
	 * get the ground-state energy
	 * @note zero-operator term in expansion of equation (20) in (Toth 2015)
	 */
	t_real CalcGroundStateEnergy() const
	{
		t_real E = 0.;

		for(const ExchangeTerm& term : GetExchangeTerms())
		{
			// check if the site indices are valid
			if(!CheckMagneticSite(term.site1_calc) || !CheckMagneticSite(term.site2_calc))
				continue;

			const MagneticSite& s_i = GetMagneticSite(term.site1_calc);
			const MagneticSite& s_j = GetMagneticSite(term.site2_calc);

			const t_mat J = CalcRealJ(term);  // Q = 0 -> no rotation needed

			const t_vec Si = s_i.spin_mag_calc * s_i.trafo_z_calc;
			const t_vec Sj = s_j.spin_mag_calc * s_j.trafo_z_calc;

			E += tl2::inner_noconj<t_vec>(Si, J * Sj).real();
		}

		return E;
	}



	/**
	 * minimise energy to found ground state
	 */
#if defined(__TLIBS2_USE_MINUIT__) && defined(__TLIBS2_MAGDYN_USE_MINUIT__)
	bool CalcGroundState(const std::unordered_set<std::string>* fixed_params = nullptr,
		bool verbose = false)
	{
		// function to minimise the state's energy
		auto func = [this](const std::vector<tl2::t_real_min>& args)
		{
			// set new spin configuration
			auto dyn = *this;

			for(t_size site_idx = 0; site_idx < dyn.GetMagneticSitesCount(); ++site_idx)
			{
				MagneticSite& site = dyn.m_sites[site_idx];

				t_real u = args[site_idx * 2 + 0];
				t_real v = args[site_idx * 2 + 1];

				const auto [ phi, theta ] = tl2::uv_to_sph<t_real>(u, v);
				const auto [ x, y, z ] = tl2::sph_to_cart<t_real>(1., phi, theta);

				site.spin_dir[0] = tl2::var_to_str(x, m_prec);
				site.spin_dir[1] = tl2::var_to_str(y, m_prec);
				site.spin_dir[2] = tl2::var_to_str(z, m_prec);

				dyn.CalcMagneticSite(site);
			}

			// ground state energy with the new configuration
			return dyn.CalcGroundStateEnergy();
		};

		// add minimisation parameters and initial values
		t_size num_args = GetMagneticSitesCount() * 2;
		std::vector<std::string> params;
		std::vector<t_real> vals, errs;
		std::vector<t_real> lower_lims, upper_lims;
		std::vector<bool> fixed;
		params.reserve(num_args);
		vals.reserve(num_args);
		errs.reserve(num_args);
		lower_lims.reserve(num_args);
		upper_lims.reserve(num_args);
		fixed.reserve(num_args);

		for(const MagneticSite& site : GetMagneticSites())
		{
			const t_vec_real& S = site.spin_dir_calc;
			const auto [ rho, phi, theta ] =  tl2::cart_to_sph<t_real>(S[0], S[1], S[2]);
			const auto [ u, v ] = tl2::sph_to_uv<t_real>(phi, theta);

			std::string phi_name = site.name + "_phi";     // phi or u parameter
			std::string theta_name = site.name + "_theta"; // theta or v parameter
			params.push_back(phi_name);
			params.push_back(theta_name);

			bool fix_phi = false;
			bool fix_theta = false;
			if(fixed_params)
			{
				fix_phi = (fixed_params->find(phi_name) != fixed_params->end());
				fix_theta = (fixed_params->find(theta_name) != fixed_params->end());
			}
			fixed.push_back(fix_phi);
			fixed.push_back(fix_theta);

			vals.push_back(u);
			vals.push_back(v);

			lower_lims.push_back(0. -m_eps);
			lower_lims.push_back(0. -m_eps);

			upper_lims.push_back(1. + m_eps);
			upper_lims.push_back(1. + m_eps);

			errs.push_back(0.1);
			errs.push_back(0.1);
		}

		if(tl2::minimise_dynargs<t_real>(num_args, func,
			params, vals, errs, &fixed, &lower_lims, &upper_lims, verbose))
		{
			// set the spins to the newly-found ground state
			for(t_size site_idx = 0; site_idx < GetMagneticSitesCount(); ++site_idx)
			{
				MagneticSite& site = m_sites[site_idx];

				t_real u = vals[site_idx * 2 + 0];
				t_real v = vals[site_idx * 2 + 1];
				tl2::set_eps_round<t_real>(u, m_eps);
				tl2::set_eps_round<t_real>(v, m_eps);

				const auto [ phi, theta ] = tl2::uv_to_sph<t_real>(u, v);
				auto [ x, y, z ] = tl2::sph_to_cart<t_real>(1., phi, theta);
				tl2::set_eps_round<t_real>(x, m_eps);
				tl2::set_eps_round<t_real>(y, m_eps);
				tl2::set_eps_round<t_real>(z, m_eps);

				site.spin_dir[0] = tl2::var_to_str(x, m_prec);
				site.spin_dir[1] = tl2::var_to_str(y, m_prec);
				site.spin_dir[2] = tl2::var_to_str(z, m_prec);

				CalcMagneticSite(site);

#ifdef __TLIBS2_MAGDYN_DEBUG_OUTPUT__
				using namespace tl2_ops;
				std::cout << site.name
					<< ": u = " << u << ", "
					<< "v = " << v << ", "
					<< "phi = " << phi << ", "
					<< "theta = " << theta << ", "
					<< "S = " << site.spin_dir_calc << std::endl;
#endif
			}

			return true;
		}
		else
		{
			std::cerr << "Magdyn error: Ground state minimisation did not converge."
				<< std::endl;
			return false;
		}
	}
#else  // __TLIBS2_USE_MINUIT__
	bool CalcGroundState(const std::unordered_set<std::string>* = nullptr, bool = false)
	{
		std::cerr << "Magdyn error: Ground state minimisation support disabled."
			<< std::endl;
		return false;
	}
#endif

	// --------------------------------------------------------------------



	// --------------------------------------------------------------------
	// loading and saving
	// --------------------------------------------------------------------
	/**
	 * generates the dispersion plot along the given q path
	 */
	void SaveDispersion(const std::string& filename,
		t_real h_start, t_real k_start, t_real l_start,
		t_real h_end, t_real k_end, t_real l_end,
		t_size num_qs = 128, t_size num_threads = 4) const
	{
		std::ofstream ofstr{filename};
		SaveDispersion(ofstr,
			h_start, k_start, l_start,
			h_end, k_end, l_end, num_qs,
			num_threads);
	}



	/**
	 * generates the dispersion plot along the given q path
	 */
	void SaveDispersion(std::ostream& ostr,
		t_real h_start, t_real k_start, t_real l_start,
		t_real h_end, t_real k_end, t_real l_end,
		t_size num_qs = 128, t_size num_threads = 4) const
	{
		ostr.precision(m_prec);

		ostr	<< std::setw(m_prec*2) << std::left << "# h"
			<< std::setw(m_prec*2) << std::left << "k"
			<< std::setw(m_prec*2) << std::left << "l"
			<< std::setw(m_prec*2) << std::left << "E"
			<< std::setw(m_prec*2) << std::left << "w"
			<< std::setw(m_prec*2) << std::left << "S_xx"
			<< std::setw(m_prec*2) << std::left << "S_yy"
			<< std::setw(m_prec*2) << std::left << "S_zz"
			<< std::endl;

		SofQEs results = CalcDispersion(h_start, k_start, l_start,
			h_end, k_end, l_end, num_qs, num_threads);

		// print results
		for(const auto& result : results)
		{
			for(const EnergyAndWeight& E_and_S : result.E_and_S)
			{
				ostr	<< std::setw(m_prec*2) << std::left << result.h
					<< std::setw(m_prec*2) << std::left << result.k
					<< std::setw(m_prec*2) << std::left << result.l
					<< std::setw(m_prec*2) << E_and_S.E
					<< std::setw(m_prec*2) << E_and_S.weight
					<< std::setw(m_prec*2) << E_and_S.S_perp(0, 0).real()
					<< std::setw(m_prec*2) << E_and_S.S_perp(1, 1).real()
					<< std::setw(m_prec*2) << E_and_S.S_perp(2, 2).real()
					<< std::endl;
			}
		}
	}



	/**
	 * load a configuration from a file
	 */
	bool Load(const std::string& filename)
	{
		try
		{
			// properties tree
			boost::property_tree::ptree node;

			// read xml file
			std::ifstream ifstr{filename};
			boost::property_tree::read_xml(ifstr, node);

			const auto &magdyn = node.get_child("magdyn");
			return Load(magdyn);
		}
		catch(const std::exception& ex)
		{
			std::cerr << "Magdyn error: Could not load \"" << filename << "\"."
				<< " Reason: " << ex.what()
				<< std::endl;

			return false;
		}
	}



	/**
	 * save the configuration to a file
	 */
	bool Save(const std::string& filename) const
	{
		// properties tree
		boost::property_tree::ptree node;

		if(!Save(node))
			return false;

		boost::property_tree::ptree root_node;
		root_node.put_child("magdyn", node);

		// write xml file
		std::ofstream ofstr{filename};
		if(!ofstr)
			return false;

		ofstr.precision(m_prec);
		boost::property_tree::write_xml(ofstr, root_node,
			boost::property_tree::xml_writer_make_settings(
				'\t', 1, std::string{"utf-8"}));
		return true;
	}



	/**
	 * load a configuration from a property tree
	 */
	bool Load(const boost::property_tree::ptree& node)
	{
		// check signature
		if(auto optInfo = node.get_optional<std::string>("meta.info");
			!optInfo || !(*optInfo == std::string{"magdyn_tool"}))
		{
			return false;
		}

		Clear();

		// variables
		if(auto vars = node.get_child_optional("variables"); vars)
		{
			m_variables.reserve(vars->size());

			for(const auto &var : *vars)
			{
				auto name = var.second.get_optional<std::string>("name");
				if(!name)
					continue;

				Variable variable;
				variable.name = *name;
				variable.value = var.second.get<t_cplx>("value", 0.);

				AddVariable(std::move(variable));
			}
		}

		// magnetic sites
		if(auto sites = node.get_child_optional("atom_sites"); sites)
		{
			m_sites.reserve(sites->size());
			std::unordered_set<std::string> seen_names;
			t_size unique_name_counter = 1;

			for(const auto &site : *sites)
			{
				MagneticSite magnetic_site;

				magnetic_site.name = site.second.get<std::string>("name", "");
				if(magnetic_site.name == "")
					magnetic_site.name = "site_" + tl2::var_to_str(GetMagneticSitesCount());

				if(seen_names.find(magnetic_site.name) != seen_names.end())
				{
					// try to create a unique name
					magnetic_site.name += "_" + tl2::var_to_str(unique_name_counter, m_prec);
					++unique_name_counter;
				}
				else
				{
					seen_names.insert(magnetic_site.name);
				}

				magnetic_site.pos_calc = tl2::create<t_vec_real>(
				{
					site.second.get<t_real>("position_x", 0.),
					site.second.get<t_real>("position_y", 0.),
					site.second.get<t_real>("position_z", 0.),
				});

				magnetic_site.sym_idx = site.second.get<t_size>("symmetry_index", 0);

				magnetic_site.pos[0] = site.second.get<std::string>("position_x", "0");
				magnetic_site.pos[1] = site.second.get<std::string>("position_y", "0");
				magnetic_site.pos[2] = site.second.get<std::string>("position_z", "0");

				magnetic_site.spin_dir[0] = site.second.get<std::string>("spin_x", "0");
				magnetic_site.spin_dir[1] = site.second.get<std::string>("spin_y", "0");
				magnetic_site.spin_dir[2] = site.second.get<std::string>("spin_z", "1");

				magnetic_site.spin_ortho[0] = site.second.get<std::string>("spin_ortho_x", "");
				magnetic_site.spin_ortho[1] = site.second.get<std::string>("spin_ortho_y", "");
				magnetic_site.spin_ortho[2] = site.second.get<std::string>("spin_ortho_z", "");

				magnetic_site.spin_mag = site.second.get<std::string>("spin_magnitude", "1");

				if(magnetic_site.g_e.size1() == 0 || magnetic_site.g_e.size2() == 0)
					magnetic_site.g_e = tl2::g_e<t_real> * tl2::unit<t_mat>(3);

				for(std::uint8_t i = 0; i < 3; ++i)
				for(std::uint8_t j = 0; j < 3; ++j)
				{
					std::string g_name = std::string{"gfactor_"} +
						std::string{g_comp_names[i]} +
						std::string{g_comp_names[j]};

					if(auto g_comp = site.second.get_optional<t_cplx>(g_name); g_comp)
						magnetic_site.g_e(i, j) = *g_comp;
				}

				AddMagneticSite(std::move(magnetic_site));
			}
		}

		// exchange terms / couplings
		if(auto terms = node.get_child_optional("exchange_terms"); terms)
		{
			m_exchange_terms.reserve(terms->size());
			std::unordered_set<std::string> seen_names;
			t_size unique_name_counter = 1;

			for(const auto &term : *terms)
			{
				ExchangeTerm exchange_term;

				exchange_term.name = term.second.get<std::string>("name", "");
				if(exchange_term.name == "")
					exchange_term.name = "coupling_" + tl2::var_to_str(GetExchangeTermsCount());

				if(seen_names.find(exchange_term.name) != seen_names.end())
				{
					// try to create a unique name
					exchange_term.name += "_" + tl2::var_to_str(unique_name_counter, m_prec);
					++unique_name_counter;
				}
				else
				{
					seen_names.insert(exchange_term.name);
				}

				exchange_term.site1_calc = term.second.get<t_size>("atom_1_index", 0);
				exchange_term.site2_calc = term.second.get<t_size>("atom_2_index", 0);

				if(auto name1 = term.second.get_optional<std::string>("atom_1_name"); name1)
				{
					// get the magnetic site index via the name
					if(const MagneticSite* sites1 = FindMagneticSite(*name1); sites1)
						exchange_term.site1 = sites1->name;
					else
					{
						std::cerr << "Magdyn error: Site 1 name \"" << *name1 << "\" "
							<< "was not found in coupling \"" << exchange_term.name << "\"."
							<< std::endl;
					}
				}
				else
				{
					// get the magnetic site name via the index
					exchange_term.site1 = GetMagneticSite(exchange_term.site1_calc).name;
				}

				if(auto name2 = term.second.get_optional<std::string>("atom_2_name"); name2)
				{
					if(const MagneticSite* sites2 = FindMagneticSite(*name2); sites2)
						exchange_term.site2 = sites2->name;
					else
					{
						std::cerr << "Magdyn error: Site 2 name \"" << *name2 << "\" "
							<< "was not found in coupling \"" << exchange_term.name << "\"."
							<< std::endl;
					}
				}
				else
				{
					// get the magnetic site name via the index
					exchange_term.site2 = GetMagneticSite(exchange_term.site2_calc).name;
				}

				exchange_term.dist_calc = tl2::create<t_vec_real>(
				{
					term.second.get<t_real>("distance_x", 0.),
					term.second.get<t_real>("distance_y", 0.),
					term.second.get<t_real>("distance_z", 0.),
				});

				exchange_term.dist[0] = term.second.get<std::string>("distance_x", "0");
				exchange_term.dist[1] = term.second.get<std::string>("distance_y", "0");
				exchange_term.dist[2] = term.second.get<std::string>("distance_z", "0");

				exchange_term.sym_idx = term.second.get<t_size>("symmetry_index", 0);

				exchange_term.J = term.second.get<std::string>("interaction", "0");

				for(std::uint8_t i = 0; i < 3; ++i)
				{
					exchange_term.dmi[i] = term.second.get<std::string>(
						std::string{"dmi_"} + std::string{g_comp_names[i]}, "0");

					for(std::uint8_t j = 0; j < 3; ++j)
					{
						exchange_term.Jgen[i][j] = term.second.get<std::string>(
							std::string{"gen_"} +
							std::string{g_comp_names[i]} +
							std::string{g_comp_names[j]}, "0");
					}
				}

				AddExchangeTerm(std::move(exchange_term));
			}
		}

		// external field
		if(auto field = node.get_child_optional("field"); field)
		{
			ExternalField thefield;

			thefield.mag = 0.;
			thefield.align_spins = false;

			thefield.dir = tl2::create<t_vec_real>(
			{
				field->get<t_real>("direction_h", 0.),
				field->get<t_real>("direction_k", 0.),
				field->get<t_real>("direction_l", 1.),
			});

			if(auto optVal = field->get_optional<t_real>("magnitude"))
				thefield.mag = *optVal;

			if(auto optVal = field->get_optional<bool>("align_spins"))
				thefield.align_spins = *optVal;

			SetExternalField(thefield);
		}

		// temperature
		m_temperature = node.get<t_real>("temperature", -1.);

		// form factor
		SetMagneticFormFactor(node.get<std::string>("magnetic_form_factor", ""));

		// ordering vector
		if(auto ordering = node.get_child_optional("ordering"); ordering)
		{
			t_vec_real ordering_vec = tl2::create<t_vec_real>(
			{
				ordering->get<t_real>("h", 0.),
				ordering->get<t_real>("k", 0.),
				ordering->get<t_real>("l", 0.),
			});

			SetOrderingWavevector(ordering_vec);
		}

		// rotation axis
		if(auto axis = node.get_child_optional("rotation_axis"); axis)
		{
			t_vec_real rotaxis = tl2::create<t_vec_real>(
			{
				axis->get<t_real>("h", 1.),
				axis->get<t_real>("k", 0.),
				axis->get<t_real>("l", 0.),
			});

			SetRotationAxis(rotaxis);
		}

		// crystal lattice
		t_real a = node.get<t_real>("xtal.a", 5.);
		t_real b = node.get<t_real>("xtal.b", 5.);
		t_real c = node.get<t_real>("xtal.c", 5.);
		t_real alpha = node.get<t_real>("xtal.alpha", 90.) / 180. * tl2::pi<t_real>;
		t_real beta = node.get<t_real>("xtal.beta", 90.) / 180. * tl2::pi<t_real>;
		t_real gamma = node.get<t_real>("xtal.gamma", 90.) / 180. * tl2::pi<t_real>;
		SetCrystalLattice(a, b, c, alpha, beta, gamma);

		// scattering plane
		t_real ah = node.get<t_real>("xtal.plane_ah", 1.);
		t_real ak = node.get<t_real>("xtal.plane_ak", 0.);
		t_real al = node.get<t_real>("xtal.plane_al", 0.);
		t_real bh = node.get<t_real>("xtal.plane_bh", 0.);
		t_real bk = node.get<t_real>("xtal.plane_bk", 1.);
		t_real bl = node.get<t_real>("xtal.plane_bl", 0.);
		SetScatteringPlane(ah, ak, al, bh, bk, bl);

		CalcExternalField();
		CalcMagneticSites();
		CalcExchangeTerms();

		return true;
	}



	/**
	 * save the configuration to a property tree
	 */
	bool Save(boost::property_tree::ptree& node) const
	{
		// write signature
		node.put<std::string>("meta.info", "magdyn_tool");
		node.put<std::string>("meta.date", tl2::epoch_to_str<t_real>(tl2::epoch<t_real>()));
		node.put<std::string>("meta.doi_tlibs", "https://doi.org/10.5281/zenodo.5717779");

		// external field
		if(m_field.dir.size() == 3)
		{
			node.put<t_real>("field.direction_h", m_field.dir[0]);
			node.put<t_real>("field.direction_k", m_field.dir[1]);
			node.put<t_real>("field.direction_l", m_field.dir[2]);
		}
		node.put<t_real>("field.magnitude", m_field.mag);
		node.put<bool>("field.align_spins", m_field.align_spins);

		// ordering vector
		if(m_ordering.size() == 3)
		{
			node.put<t_real>("ordering.h", m_ordering[0]);
			node.put<t_real>("ordering.k", m_ordering[1]);
			node.put<t_real>("ordering.l", m_ordering[2]);
		}

		// rotation axis
		if(m_rotaxis.size() == 3)
		{
			node.put<t_real>("rotation_axis.h", m_rotaxis[0]);
			node.put<t_real>("rotation_axis.k", m_rotaxis[1]);
			node.put<t_real>("rotation_axis.l", m_rotaxis[2]);
		}

		// temperature
		node.put<t_real>("temperature", m_temperature);

		// form factor
		node.put<std::string>("magnetic_form_factor", GetMagneticFormFactor());

		// variables
		for(const Variable& var : GetVariables())
		{
			boost::property_tree::ptree itemNode;
			itemNode.put<std::string>("name", var.name);
			itemNode.put<t_cplx>("value", var.value);

			node.add_child("variables.variable", itemNode);
		}

		// magnetic sites
		for(const MagneticSite& site : GetMagneticSites())
		{
			boost::property_tree::ptree itemNode;
			itemNode.put<std::string>("name", site.name);

			itemNode.put<std::string>("position_x", site.pos[0]);
			itemNode.put<std::string>("position_y", site.pos[1]);
			itemNode.put<std::string>("position_z", site.pos[2]);

			itemNode.put<t_size>("symmetry_index", site.sym_idx);

			itemNode.put<std::string>("spin_x", site.spin_dir[0]);
			itemNode.put<std::string>("spin_y", site.spin_dir[1]);
			itemNode.put<std::string>("spin_z", site.spin_dir[2]);

			itemNode.put<std::string>("spin_ortho_x", site.spin_ortho[0]);
			itemNode.put<std::string>("spin_ortho_y", site.spin_ortho[1]);
			itemNode.put<std::string>("spin_ortho_z", site.spin_ortho[2]);

			itemNode.put<std::string>("spin_magnitude", site.spin_mag);

			for(std::uint8_t i = 0; i < site.g_e.size1(); ++i)
			for(std::uint8_t j = 0; j < site.g_e.size2(); ++j)
			{
				itemNode.put<t_cplx>(std::string{"gfactor_"} +
					std::string{g_comp_names[i]} +
					std::string{g_comp_names[j]}, site.g_e(i, j));
			}

			node.add_child("atom_sites.site", itemNode);
		}

		// exchange terms
		for(const ExchangeTerm& term : GetExchangeTerms())
		{
			boost::property_tree::ptree itemNode;
			itemNode.put<std::string>("name", term.name);

			// save the magnetic site names and indices
			itemNode.put<t_size>("atom_1_index", term.site1_calc);
			itemNode.put<t_size>("atom_2_index", term.site2_calc);
			itemNode.put<std::string>("atom_1_name", term.site1);
			itemNode.put<std::string>("atom_2_name", term.site2);

			itemNode.put<std::string>("distance_x", term.dist[0]);
			itemNode.put<std::string>("distance_y", term.dist[1]);
			itemNode.put<std::string>("distance_z", term.dist[2]);

			itemNode.put<t_size>("symmetry_index", term.sym_idx);

			itemNode.put<std::string>("interaction", term.J);

			for(std::uint8_t i = 0; i < 3; ++i)
			{
				itemNode.put<std::string>(std::string{"dmi_"} +
					std::string{g_comp_names[i]}, term.dmi[i]);

				for(std::uint8_t j = 0; j < 3; ++j)
				{
					itemNode.put<std::string>(std::string{"gen_"} +
						std::string{g_comp_names[i]} +
						std::string{g_comp_names[j]},
						term.Jgen[i][j]);
				}
			}

			node.add_child("exchange_terms.term", itemNode);
		}

		// crystal lattice
		node.put<t_real>("xtal.a", m_xtallattice[0]);
		node.put<t_real>("xtal.b", m_xtallattice[1]);
		node.put<t_real>("xtal.c", m_xtallattice[2]);
		node.put<t_real>("xtal.alpha", m_xtalangles[0] / tl2::pi<t_real> * 180.);
		node.put<t_real>("xtal.beta", m_xtalangles[1] / tl2::pi<t_real> * 180.);
		node.put<t_real>("xtal.gamma", m_xtalangles[2] / tl2::pi<t_real> * 180.);

		// scattering plane
		// x vector
		node.put<t_real>("xtal.plane_ah", m_scatteringplane[0][0]);
		node.put<t_real>("xtal.plane_ak", m_scatteringplane[0][1]);
		node.put<t_real>("xtal.plane_al", m_scatteringplane[0][2]);
		// y vector
		node.put<t_real>("xtal.plane_bh", m_scatteringplane[1][0]);
		node.put<t_real>("xtal.plane_bk", m_scatteringplane[1][1]);
		node.put<t_real>("xtal.plane_bl", m_scatteringplane[1][2]);
		// up vector (saving is optional)
		node.put<t_real>("xtal.plane_ch", m_scatteringplane[2][0]);
		node.put<t_real>("xtal.plane_ck", m_scatteringplane[2][1]);
		node.put<t_real>("xtal.plane_cl", m_scatteringplane[2][2]);

		return true;
	}
	// --------------------------------------------------------------------



protected:
	/**
	 * converts the rotation matrix rotating the local spins to ferromagnetic
	 * [001] directions into the vectors comprised of the matrix columns
	 * @see equation (9) and (51) from (Toth 2015)
	 */
	std::tuple<t_vec, t_vec> rot_to_trafo(const t_mat& R)
	{
		const t_vec xy_plane = tl2::col<t_mat, t_vec>(R, 0)
			 + s_imag * tl2::col<t_mat, t_vec>(R, 1);
		const t_vec z = tl2::col<t_mat, t_vec>(R, 2);

		return std::make_tuple(xy_plane, z);
	}



	/**
	 * rotate local spin to ferromagnetic [001] direction
	 * @see equations (7) and (9) from (Toth 2015)
	 */
	std::tuple<t_vec, t_vec> spin_to_trafo(const t_vec_real& spin_dir)
	{
		const t_mat_real _rot = tl2::rotation<t_mat_real, t_vec_real>(
			spin_dir, m_zdir, &m_rotaxis, m_eps);

		const t_mat rot = tl2::convert<t_mat, t_mat_real>(_rot);
		return rot_to_trafo(rot);
	}



private:
	// magnetic sites
	MagneticSites m_sites{};

	// magnetic couplings
	ExchangeTerms m_exchange_terms{};

	// open variables in expressions
	Variables m_variables{};

	// external field
	ExternalField m_field{};
	// matrix to rotate field into the [001] direction
	t_mat m_rot_field{ tl2::unit<t_mat>(3) };

	// ordering wave vector for incommensurate structures
	t_vec_real m_ordering{ tl2::zero<t_vec_real>(3) };

	// helix rotation axis for incommensurate structures
	t_vec_real m_rotaxis{ tl2::create<t_vec_real>({ 1., 0., 0. }) };

	// calculate the hamiltonian for Q, Q+ordering, and Q-ordering
	bool m_calc_H{ true };
	bool m_calc_Hp{ true };
	bool m_calc_Hm{ true };

	// direction to rotation spins into, usually [001]
	t_vec_real m_zdir{ tl2::create<t_vec_real>({ 0., 0., 1. }) };

	// temperature (-1: disable bose factor)
	t_real m_temperature{ -1. };

	// bose cutoff energy to avoid infinities
	t_real m_bose_cutoff{ 0.025 };

	// formula for the magnetic form factor
	std::string m_magffact_formula{};
	tl2::ExprParser<t_cplx> m_magffact{};

	// crystal lattice
	t_real m_xtallattice[3]{ 5., 5., 5. };
	t_real m_xtalangles[3]
	{
		t_real(0.5) * tl2::pi<t_real>,
		t_real(0.5) * tl2::pi<t_real>,
		t_real(0.5) * tl2::pi<t_real>
	};
	t_mat_real m_xtalA{ tl2::unit<t_mat_real>(3) };
	t_mat_real m_xtalB{ tl2::unit<t_mat_real>(3) };
	t_mat_real m_xtalUB{ tl2::unit<t_mat_real>(3) };
	t_mat_real m_xtalUBinv{ tl2::unit<t_mat_real>(3) };

	//scattering plane
	t_vec_real m_scatteringplane[3]
	{
		tl2::create<t_vec_real>({ 1., 0., 0. }),  // in-plane, x
		tl2::create<t_vec_real>({ 0., 1., 0. }),  // in-plane, y
		tl2::create<t_vec_real>({ 0., 0., 1. }),  // out-of-plane, z
	};

	// settings
	bool m_is_incommensurate{ false };
	bool m_force_incommensurate{ false };
	bool m_unite_degenerate_energies{ true };
	bool m_perform_checks{ true };

	// settings for cholesky decomposition
	t_size m_tries_chol{ 50 };
	t_real m_delta_chol{ 0.0025 };

	// precisions
	t_real m_eps{ 1e-6 };
	int m_prec{ 6 };

	// conventions
	t_real m_phase_sign{ -1. };

	// constants
	static constexpr const t_cplx s_imag { t_real(0), t_real(1) };
	static constexpr const t_real s_twopi { t_real(2)*tl2::pi<t_real> };
	static constexpr const std::array<std::string_view, 3> g_comp_names{{ "x", "y", "z" }};
};

}
#endif
