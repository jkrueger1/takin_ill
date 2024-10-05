/**
 * tlibs2 -- fitting and minimisation library
 * @author Tobias Weber <tobias.weber@tum.de>, <tweber@ill.fr>
 * @date 2012-2024
 * @note Forked on 7-Nov-2018 from my privately and TUM-PhD-developed "tlibs" project (https://github.com/t-weber/tlibs).
 * @license GPLv3, see 'LICENSE' file
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

#ifndef __TLIBS2_FITTER_H__
#define __TLIBS2_FITTER_H__

#if __has_include(<Minuit2/FCNBase.h>) && __has_include(<Minuit2/MnTraceObject.h>)
	#include <Minuit2/FCNBase.h>
	#include <Minuit2/MnFcn.h>
	#include <Minuit2/FunctionMinimum.h>
	#include <Minuit2/MnMigrad.h>
	#include <Minuit2/MnPrint.h>

	#define __TLIBS2_USE_MINUIT__
#else
	//#pragma message("tlibs2: Disabling Minuit library (not found).")
#endif

#include <vector>
#include <iostream>
#include <string>
#include <algorithm>
#include <type_traits>
#include <exception>

#include "expr.h"
#include "maths.h"


namespace tl2 {


// ----------------------------------------------------------------------------
// stop request handling
// ----------------------------------------------------------------------------
struct StopRequestException : public std::runtime_error
{
	StopRequestException(const char* msg = "") : runtime_error{msg}
	{}
};


class StopRequest
{
public:
	StopRequest() = default;
	virtual ~StopRequest() = default;


	StopRequest(const StopRequest& req)
		: m_stop_requested{req.m_stop_requested}
	{
	}


	const StopRequest& operator=(const StopRequest& req)
	{
		m_stop_requested = req.m_stop_requested;
		return *this;
	}


	void SetStopRequest(bool *b)
	{
		m_stop_requested = b;
	}


	void HandleStopRequest() const
	{
		if(!m_stop_requested || !*m_stop_requested)
			return;

		// the only way to get out of an ongoing minuit operation
		throw StopRequestException("Stop requested.");
	}


private:
	bool *m_stop_requested{};
};
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// Minuit interface
// @see http://seal.cern.ch/documents/minuit/mnusersguide.pdf
// ----------------------------------------------------------------------------

#ifdef __TLIBS2_USE_MINUIT__
using t_real_min = std::invoke_result_t<
	decltype(&ROOT::Minuit2::MnFcn::Up), ROOT::Minuit2::MnFcn>;



template<class t_real>
class FitterFuncModel
{
public:
	virtual ~FitterFuncModel() = default;

	virtual bool SetParams(const std::vector<t_real>& vecParams) = 0;
	virtual t_real operator()(t_real x) const = 0;
	virtual FitterFuncModel<t_real>* copy() const = 0;
};



/**
 * interface using supplied functions with a fixed number of parameters
 * num_args also includes the "x" parameter to the function, m_vecVals does not
 */
template<class t_real, std::size_t num_args, typename t_func>
class FitterLamFuncModel : public FitterFuncModel<t_real>
{
protected:
	t_func m_func{};
	std::vector<t_real> m_vecVals{};
	bool m_bSeparateFreeParam{true};  // separate "x" from parameters (for fitter)

public:
	FitterLamFuncModel(t_func func, bool bSeparateX = true)
		: m_func{func}, m_vecVals{}, m_bSeparateFreeParam{bSeparateX}
	{
		m_vecVals.resize(m_bSeparateFreeParam ? num_args - 1 : num_args);
	}


	virtual bool SetParams(const std::vector<t_real>& vecParams) override
	{
		for(std::size_t i = 0; i < std::min(vecParams.size(), m_vecVals.size()); ++i)
			m_vecVals[i] = vecParams[i];
		return true;
	}


	virtual t_real operator()(t_real x = t_real(0)) const override
	{
		std::vector<t_real> vecValsWithX;
		vecValsWithX.reserve(num_args);

		if(m_bSeparateFreeParam)
		{
			vecValsWithX.push_back(x);
			for(t_real d : m_vecVals)
				vecValsWithX.push_back(d);
		}

		const std::vector<t_real> *pvecVals = m_bSeparateFreeParam ? &vecValsWithX : &m_vecVals;
		t_real funcval = call<num_args, t_func, t_real, std::vector>(m_func, *pvecVals);
		return funcval;
	}


	virtual FitterLamFuncModel* copy() const override
	{
		FitterLamFuncModel<t_real, num_args, t_func>* pMod =
			new FitterLamFuncModel<t_real, num_args, t_func>(m_func);

		pMod->m_vecVals = this->m_vecVals;
		pMod->m_bSeparateFreeParam = this->m_bSeparateFreeParam;

		return pMod;
	}
};



/**
 * interface using supplied functions with a dynamic number of parameters
 */
template<class t_real, typename t_func>
class FitterDynLamFuncModel : public FitterFuncModel<t_real>
{
protected:
	std::size_t m_num_args{1};
	t_func m_func{};
	std::vector<t_real> m_vecVals{};
	bool m_bSeparateFreeParam{true};  // separate "x" from parameters (for fitter)

public:
	/**
	 * num_args also includes the "x" parameter to the function, m_vecVals does not
	 */
	FitterDynLamFuncModel(std::size_t num_args, t_func func, bool bSeparateX = true)
	: m_num_args{num_args}, m_func{func}, m_vecVals{}, m_bSeparateFreeParam{bSeparateX}
	{
		m_vecVals.resize(m_bSeparateFreeParam ? m_num_args - 1 : m_num_args);
	}


	virtual bool SetParams(const std::vector<t_real>& vecParams) override
	{
		for(std::size_t i = 0; i < std::min(vecParams.size(), m_vecVals.size()); ++i)
			m_vecVals[i] = vecParams[i];
		return true;
	}


	virtual t_real operator()(t_real x = t_real(0)) const override
	{
		std::vector<t_real> vecValsWithX;
		vecValsWithX.reserve(m_num_args);

		if(m_bSeparateFreeParam)
		{
			vecValsWithX.push_back(x);
			for(t_real d : m_vecVals)
				vecValsWithX.push_back(d);
		}

		return m_func(m_bSeparateFreeParam ? vecValsWithX : m_vecVals);
	}


	virtual FitterDynLamFuncModel* copy() const override
	{
		FitterDynLamFuncModel<t_real, t_func>* pMod =
		new FitterDynLamFuncModel<t_real, t_func>(m_num_args, m_func);

		pMod->m_vecVals = this->m_vecVals;
		pMod->m_bSeparateFreeParam = this->m_bSeparateFreeParam;

		return pMod;
	}
};



/**
 * interface using supplied functions
 * num_args also includes the "x" parameter to the function, m_vecVals does not
 */
template<class t_real>
class FitterParsedFuncModel : public FitterFuncModel<t_real>
{
protected:
	std::string m_func;

	std::string m_xName = "x";
	const std::vector<std::string>& m_vecNames;
	std::vector<t_real> m_vecVals;

	ExprParser<t_real> m_expr{};


public:
	FitterParsedFuncModel(const std::string& func, const std::string& xName,
		const std::vector<std::string>& vecNames)
		: m_func{func}, m_xName{xName}, m_vecNames{vecNames}
	{
		if(!m_expr.parse(m_func))
			throw std::runtime_error("Could not parse function.");
	}


	virtual bool SetParams(const std::vector<t_real>& vecParams) override
	{
		m_vecVals.resize(vecParams.size());
		for(std::size_t i = 0; i < std::min(vecParams.size(), m_vecVals.size()); ++i)
			m_vecVals[i] = vecParams[i];
		return true;
	}


	virtual t_real operator()(t_real x = t_real(0)) const override
	{
		// copy the parsed expression to be thread safe
		ExprParser<t_real> expr = m_expr;

		// x is not used for minimiser
		if(m_xName != "")
			expr.register_var(m_xName, x);

		for(std::size_t i = 0; i < m_vecVals.size(); ++i)
			expr.register_var(m_vecNames[i], m_vecVals[i]);

		t_real val = expr.eval();
		return val;
	}


	virtual FitterParsedFuncModel* copy() const override
	{
		return new FitterParsedFuncModel<t_real>(m_func, m_xName, m_vecNames);
	}
};

// ----------------------------------------------------------------------------



/**
 * generic chi^2 calculation for fitting
 * @see http://seal.cern.ch/documents/minuit/mnusersguide.pdf
 */
template<class t_real = t_real_min>
class Chi2Function : public ROOT::Minuit2::FCNBase, public StopRequest
{
protected:
	const FitterFuncModel<t_real_min> *m_pfkt = nullptr;

	std::size_t m_num_pts = 0;
	const t_real* m_px = nullptr;
	const t_real* m_py = nullptr;
	const t_real* m_pdy = nullptr;

	t_real_min m_dSigma = 1.;
	bool m_bDebug = false;


public:
	Chi2Function(const FitterFuncModel<t_real_min> *fkt = nullptr,
		std::size_t num_pts = 0, const t_real *px = nullptr,
		const t_real *py = nullptr, const t_real *pdy = nullptr)
		: m_pfkt{fkt}, m_num_pts{num_pts}, m_px{px}, m_py{py}, m_pdy{pdy}
	{}

	virtual ~Chi2Function() = default;


	const Chi2Function<t_real>& operator=(const Chi2Function<t_real>& other)
	{
		StopRequest::operator=(*this);

		this->m_pfkt = other.m_pfkt;
		this->m_px = other.m_px;
		this->m_py = other.m_py;
		this->m_pdy = other.m_pdy;
		this->m_dSigma = other.m_dSigma;
		this->m_bDebug = other.m_bDebug;

		return *this;
	}

	Chi2Function(const Chi2Function<t_real>& other)
	{
		operator=(other);
	}


	/*
	 * chi^2 calculation
	 * based on the example in the Minuit user's guide:
	 * http://seal.cern.ch/documents/minuit/mnusersguide.pdf
	 */
	t_real_min chi2(const std::vector<t_real_min>& vecParams) const
	{
		// cannot operate on m_pfkt directly, because Minuit
		// uses more than one thread!
		std::unique_ptr<FitterFuncModel<t_real_min>> uptrFkt(m_pfkt->copy());
		FitterFuncModel<t_real_min>* pfkt = uptrFkt.get();

		pfkt->SetParams(vecParams);
		return tl2::chi2<t_real_min, decltype(*pfkt), const t_real*>(
			*pfkt, m_num_pts, m_px, m_py, m_pdy);
	}

	virtual t_real_min Up() const override
	{
		return m_dSigma*m_dSigma;
	}

	virtual t_real_min operator()(const std::vector<t_real_min>& vecParams) const override
	{
		HandleStopRequest();

		t_real_min dChi2 = chi2(vecParams);
		if(m_bDebug)
			std::cerr << "Fitter: chi2 = " << dChi2 << "." << std::endl;
		return dChi2;
	}

	void SetSigma(t_real_min dSig)
	{
		m_dSigma = dSig;
	}

	t_real_min GetSigma() const
	{
		return m_dSigma;
	}

	void SetDebug(bool b)
	{
		m_bDebug = b;
	}
};



/**
 * function adaptor for minimisation
 * @see http://seal.cern.ch/documents/minuit/mnusersguide.pdf
 */
template<class t_real = t_real_min>
class MiniFunction : public ROOT::Minuit2::FCNBase, public StopRequest
{
protected:
	const FitterFuncModel<t_real_min> *m_pfkt = nullptr;
	t_real_min m_dSigma = 1.;

public:
	MiniFunction(const FitterFuncModel<t_real_min>* fkt = nullptr)
		: m_pfkt(fkt)
	{}

	virtual ~MiniFunction() = default;

	MiniFunction(const MiniFunction<t_real>& other)
		: m_pfkt(other.m_pfkt), m_dSigma(other.m_dSigma)
	{}

	const MiniFunction<t_real>& operator=(const MiniFunction<t_real>& other)
	{
		StopRequest::operator=(*this);

		this->m_pfkt = other.m_pfkt;
		this->m_dSigma = other.m_dSigma;

		return *this;
	}

	virtual t_real_min Up() const override
	{
		return m_dSigma*m_dSigma;
	}

	virtual t_real_min operator()(const std::vector<t_real_min>& vecParams) const override
	{
		HandleStopRequest();

		// cannot operate on m_pfkt directly, because Minuit
		// uses more than one thread!
		std::unique_ptr<FitterFuncModel<t_real_min>> uptrFkt(m_pfkt->copy());
		FitterFuncModel<t_real_min>* pfkt = uptrFkt.get();

		pfkt->SetParams(vecParams);
		return (*pfkt)(t_real_min(0));	// "0" is an ignored dummy value here
	}

	void SetSigma(t_real_min dSig)
	{
		m_dSigma = dSig;
	}

	t_real_min GetSigma() const
	{
		return m_dSigma;
	}
};



// ----------------------------------------------------------------------------



/**
 * fit function to x,y,dy data points
 */
template<class t_real = t_real_min, std::size_t num_args, typename t_func>
bool fit(t_func&& func,

	const std::vector<t_real>& vecX,
	const std::vector<t_real>& vecY,
	const std::vector<t_real>& vecYErr,

	const std::vector<std::string>& vecParamNames,	// size: num_args-1
	std::vector<t_real>& vecVals,
	std::vector<t_real>& vecErrs,
	const std::vector<bool>* pVecFixed = nullptr,

	bool bDebug = true, bool *stop_request = nullptr)
{
	try
	{
		if(!vecX.size() || !vecY.size() || !vecYErr.size())
		{
			std::cerr << "Fitter: No data given." << std::endl;
			return false;
		}

		// check if all params are fixed
		if(pVecFixed && std::all_of(pVecFixed->begin(), pVecFixed->end(),
			[](bool b) -> bool { return b; }))
			{
				std::cerr << "Fitter: All parameters are fixed." << std::endl;
				return false;
			}

		// convert vectors if value types don't match with minuit's type
		std::vector<t_real_min> vecXConverted, vecYConverted, vecYErrConverted;
		if constexpr(!std::is_same_v<t_real, t_real_min>)
		{
			vecXConverted.reserve(vecX.size());
			vecYConverted.reserve(vecY.size());
			vecYErrConverted.reserve(vecYErr.size());

			for(t_real d : vecX)
				vecXConverted.push_back(static_cast<t_real_min>(d));
			for(t_real d : vecY)
				vecYConverted.push_back(static_cast<t_real_min>(d));
			for(t_real d : vecYErr)
				vecYErrConverted.push_back(static_cast<t_real_min>(d));
		}

		FitterLamFuncModel<t_real_min, num_args, t_func> mod(func);

		std::unique_ptr<Chi2Function<t_real_min>> chi2;
		if constexpr(std::is_same_v<t_real, t_real_min>)
		{
			chi2 = std::make_unique<Chi2Function<t_real_min>>(
				&mod, vecX.size(), vecX.data(), vecY.data(), vecYErr.data());
		}
		else if constexpr(!std::is_same_v<t_real, t_real_min>)
		{
			chi2 = std::make_unique<Chi2Function<t_real_min>>(
				&mod, vecXConverted.size(), vecXConverted.data(),
				vecYConverted.data(), vecYErrConverted.data());
		}
		chi2->SetStopRequest(stop_request);

		ROOT::Minuit2::MnUserParameters params;
		for(std::size_t param_idx = 0; param_idx < vecParamNames.size(); ++param_idx)
		{
			params.Add(vecParamNames[param_idx],
				static_cast<t_real_min>(vecVals[param_idx]),
				static_cast<t_real_min>(vecErrs[param_idx]));
			if(pVecFixed && (*pVecFixed)[param_idx])
				params.Fix(vecParamNames[param_idx]);
		}

		ROOT::Minuit2::MnMigrad migrad(*chi2, params, 2);
		ROOT::Minuit2::FunctionMinimum mini = migrad();
		bool bValidFit = mini.IsValid() && mini.HasValidParameters() && mini.UserState().IsValid();

		for(std::size_t param_idx = 0; param_idx < vecParamNames.size(); ++param_idx)
		{
			vecVals[param_idx] = static_cast<t_real>(
				mini.UserState().Value(vecParamNames[param_idx]));
			vecErrs[param_idx] = static_cast<t_real>(
				std::fabs(mini.UserState().Error(vecParamNames[param_idx])));
		}

		if(bDebug)
			std::cerr << mini << std::endl;

		return bValidFit;
	}
	catch(const tl2::StopRequestException&)
	{
		throw;
	}
	catch(const std::exception& ex)
	{
		std::cerr << "Fitter: " << ex.what() << std::endl;
	}

	return false;
}



/**
 * fit expression to x,y,dy data points
 */
template<class t_real = t_real_min>
bool fit_expr(const std::string& func,

	const std::vector<t_real>& vecX,
	const std::vector<t_real>& vecY,
	const std::vector<t_real>& vecYErr,

	const std::string& strXName,
	const std::vector<std::string>& vecParamNames,	// size: num_args-1
	std::vector<t_real>& vecVals,
	std::vector<t_real>& vecErrs,
	const std::vector<bool>* pVecFixed = nullptr,

	bool bDebug = true, bool *stop_request = nullptr)
{
	try
	{
		if(!vecX.size() || !vecY.size() || !vecYErr.size())
		{
			std::cerr << "Fitter: No data given." << std::endl;
			return false;
		}

		// check if all params are fixed
		if(pVecFixed && std::all_of(pVecFixed->begin(), pVecFixed->end(),
			[](bool b) -> bool { return b; }))
			{
				std::cerr << "Fitter: All parameters are fixed." << std::endl;
				return false;
			}

		// convert vectors if value types don't match with minuit's type
		std::vector<t_real_min> vecXConverted, vecYConverted, vecYErrConverted;
		if constexpr(!std::is_same_v<t_real, t_real_min>)
		{
			vecXConverted.reserve(vecX.size());
			vecYConverted.reserve(vecY.size());
			vecYErrConverted.reserve(vecYErr.size());

			for(t_real d : vecX)
				vecXConverted.push_back(static_cast<t_real_min>(d));
			for(t_real d : vecY)
				vecYConverted.push_back(static_cast<t_real_min>(d));
			for(t_real d : vecYErr)
				vecYErrConverted.push_back(static_cast<t_real_min>(d));
		}

		FitterParsedFuncModel<t_real_min> mod(func, strXName, vecParamNames);

		std::unique_ptr<Chi2Function<t_real_min>> chi2;
		if constexpr(std::is_same_v<t_real, t_real_min>)
		{
			chi2 = std::make_unique<Chi2Function<t_real_min>>(
				&mod, vecX.size(), vecX.data(), vecY.data(), vecYErr.data());
		}
		else if constexpr(!std::is_same_v<t_real, t_real_min>)
		{
			chi2 = std::make_unique<Chi2Function<t_real_min>>(
				&mod, vecXConverted.size(), vecXConverted.data(),
				vecYConverted.data(), vecYErrConverted.data());
		}
		chi2->SetStopRequest(stop_request);

		ROOT::Minuit2::MnUserParameters params;
		for(std::size_t param_idx = 0; param_idx < vecParamNames.size(); ++param_idx)
		{
			params.Add(vecParamNames[param_idx], static_cast<t_real_min>(vecVals[param_idx]), static_cast<t_real_min>(vecErrs[param_idx]));
			if(pVecFixed && (*pVecFixed)[param_idx])
				params.Fix(vecParamNames[param_idx]);
		}

		ROOT::Minuit2::MnMigrad migrad(*chi2, params, 2);
		ROOT::Minuit2::FunctionMinimum mini = migrad();
		bool bValidFit = mini.IsValid() && mini.HasValidParameters() && mini.UserState().IsValid();

		for(std::size_t param_idx = 0; param_idx < vecParamNames.size(); ++param_idx)
		{
			vecVals[param_idx] = static_cast<t_real>(
				mini.UserState().Value(vecParamNames[param_idx]));
			vecErrs[param_idx] = static_cast<t_real>(
				std::fabs(mini.UserState().Error(vecParamNames[param_idx])));
		}

		if(bDebug)
			std::cerr << mini << std::endl;

		return bValidFit;
	}
	catch(const tl2::StopRequestException&)
	{
		throw;
	}
	catch(const std::exception& ex)
	{
		std::cerr << "Fitter: " << ex.what() << std::endl;
	}

	return false;
}



/**
 * find function minimum using a lambda function with fixed args
 */
template<class t_real = t_real_min, std::size_t num_args, typename t_func>
bool minimise(t_func&& func, const std::vector<std::string>& vecParamNames,
	std::vector<t_real>& vecVals, std::vector<t_real>& vecErrs,
	const std::vector<bool>* pVecFixed = nullptr,
	const std::vector<t_real>* pVecLowerLimits = nullptr,
	const std::vector<t_real>* pVecUpperLimits = nullptr,
	bool bDebug = true, bool *stop_request = nullptr)
{
	try
	{
		// check if all params are fixed
		if(pVecFixed && std::all_of(pVecFixed->begin(), pVecFixed->end(),
			[](bool b) -> bool { return b; }))
			{
				std::cerr << "Fitter: All parameters are fixed." << std::endl;
				return false;
			}

		FitterLamFuncModel<t_real_min, num_args, t_func> mod(func, false);
		MiniFunction<t_real_min> minfunc(&mod);
		minfunc.SetStopRequest(stop_request);

		ROOT::Minuit2::MnUserParameters params;
		for(std::size_t param_idx = 0; param_idx < vecParamNames.size(); ++param_idx)
		{
			params.Add(vecParamNames[param_idx],
				static_cast<t_real_min>(vecVals[param_idx]),
				static_cast<t_real_min>(vecErrs[param_idx]));
			if(pVecLowerLimits && pVecUpperLimits)
				params.SetLimits(vecParamNames[param_idx], (*pVecLowerLimits)[param_idx], (*pVecUpperLimits)[param_idx]);
			else if(pVecLowerLimits && !pVecUpperLimits)
				params.SetLowerLimit(vecParamNames[param_idx], (*pVecLowerLimits)[param_idx]);
			else if(pVecUpperLimits && !pVecLowerLimits)
				params.SetUpperLimit(vecParamNames[param_idx], (*pVecUpperLimits)[param_idx]);
			if(pVecFixed && (*pVecFixed)[param_idx])
				params.Fix(vecParamNames[param_idx]);
		}

		ROOT::Minuit2::MnMigrad migrad(minfunc, params, 2);
		ROOT::Minuit2::FunctionMinimum mini = migrad();
		bool bMinimumValid = mini.IsValid() && mini.HasValidParameters() && mini.UserState().IsValid();

		for(std::size_t param_idx = 0; param_idx < vecParamNames.size(); ++param_idx)
		{
			vecVals[param_idx] = static_cast<t_real>(
				mini.UserState().Value(vecParamNames[param_idx]));
			vecErrs[param_idx] = static_cast<t_real>(
				std::fabs(mini.UserState().Error(vecParamNames[param_idx])));
		}

		if(bDebug)
			std::cerr << mini << std::endl;

		return bMinimumValid;
	}
	catch(const tl2::StopRequestException&)
	{
		throw;
	}
	catch(const std::exception& ex)
	{
		std::cerr << "Fitter: " << ex.what() << std::endl;
	}

	return false;
}



/**
 * find function minimum using a lambda function with variable args
 */
template<class t_real = t_real_min, typename t_func>
bool minimise_dynargs(std::size_t num_args, t_func&& func,
	const std::vector<std::string>& vecParamNames,
	std::vector<t_real>& vecVals, std::vector<t_real>& vecErrs,
	const std::vector<bool>* pVecFixed = nullptr,
	const std::vector<t_real>* pVecLowerLimits = nullptr,
	const std::vector<t_real>* pVecUpperLimits = nullptr,
	bool bDebug = true, bool *stop_request = nullptr)
{
	try
	{
		// check if all params are fixed
		if(pVecFixed && std::all_of(pVecFixed->begin(), pVecFixed->end(),
			[](bool b) -> bool { return b; }))
			{
				std::cerr << "Fitter: All parameters are fixed." << std::endl;
				return false;
			}

		FitterDynLamFuncModel<t_real_min, t_func> mod(num_args, func, false);
		MiniFunction<t_real_min> minfunc(&mod);
		minfunc.SetStopRequest(stop_request);

		ROOT::Minuit2::MnUserParameters params;
		for(std::size_t param_idx = 0; param_idx < vecParamNames.size(); ++param_idx)
		{
			params.Add(vecParamNames[param_idx],
				static_cast<t_real_min>(vecVals[param_idx]),
				static_cast<t_real_min>(vecErrs[param_idx]));
			if(pVecLowerLimits && pVecUpperLimits)
				params.SetLimits(vecParamNames[param_idx], (*pVecLowerLimits)[param_idx], (*pVecUpperLimits)[param_idx]);
			else if(pVecLowerLimits && !pVecUpperLimits)
				params.SetLowerLimit(vecParamNames[param_idx], (*pVecLowerLimits)[param_idx]);
			else if(pVecUpperLimits && !pVecLowerLimits)
				params.SetUpperLimit(vecParamNames[param_idx], (*pVecUpperLimits)[param_idx]);
			if(pVecFixed && (*pVecFixed)[param_idx])
				params.Fix(vecParamNames[param_idx]);
		}

		ROOT::Minuit2::MnMigrad migrad(minfunc, params, 2);
		ROOT::Minuit2::FunctionMinimum mini = migrad();
		bool bMinimumValid = mini.IsValid() && mini.HasValidParameters() && mini.UserState().IsValid();

		for(std::size_t param_idx = 0; param_idx < vecParamNames.size(); ++param_idx)
		{
			vecVals[param_idx] = static_cast<t_real>(
				mini.UserState().Value(vecParamNames[param_idx]));
			vecErrs[param_idx] = static_cast<t_real>(
				std::fabs(mini.UserState().Error(vecParamNames[param_idx])));
		}

		if(bDebug)
			std::cerr << mini << std::endl;

		return bMinimumValid;
	}
	catch(const tl2::StopRequestException&)
	{
		throw;
	}
	catch(const std::exception& ex)
	{
		std::cerr << "Fitter: " << ex.what() << std::endl;
	}

	return false;
}



/**
 * find function minimum for an expression
 */
template<class t_real = t_real_min>
bool minimise_expr(const std::string& func, const std::vector<std::string>& vecParamNames,
	std::vector<t_real>& vecVals, std::vector<t_real>& vecErrs,
	const std::vector<bool>* pVecFixed = nullptr,
	bool bDebug = true, bool *stop_request = nullptr)
{
	try
	{
		// check if all params are fixed
		if(pVecFixed && std::all_of(pVecFixed->begin(), pVecFixed->end(),
			[](bool b) -> bool { return b; }))
			{
				std::cerr << "Fitter: All parameters are fixed." << std::endl;
				return false;
			}

		FitterParsedFuncModel<t_real_min> mod(func, "", vecParamNames);
		MiniFunction<t_real_min> minfunc(&mod);
		minfunc.SetStopRequest(stop_request);

		ROOT::Minuit2::MnUserParameters params;
		for(std::size_t param_idx = 0; param_idx < vecParamNames.size(); ++param_idx)
		{
			params.Add(vecParamNames[param_idx],
				static_cast<t_real_min>(vecVals[param_idx]),
				static_cast<t_real_min>(vecErrs[param_idx]));
			if(pVecFixed && (*pVecFixed)[param_idx])
				params.Fix(vecParamNames[param_idx]);
		}

		ROOT::Minuit2::MnMigrad migrad(minfunc, params, 2);
		ROOT::Minuit2::FunctionMinimum mini = migrad();
		bool bMinimumValid = mini.IsValid() && mini.HasValidParameters() && mini.UserState().IsValid();

		for(std::size_t param_idx = 0; param_idx < vecParamNames.size(); ++param_idx)
		{
			vecVals[param_idx] = static_cast<t_real>(
				mini.UserState().Value(vecParamNames[param_idx]));
			vecErrs[param_idx] = static_cast<t_real>(
				std::fabs(mini.UserState().Error(vecParamNames[param_idx])));
		}

		if(bDebug)
			std::cerr << mini << std::endl;

		return bMinimumValid;
	}
	catch(const tl2::StopRequestException&)
	{
		throw;
	}
	catch(const std::exception& ex)
	{
		std::cerr << "Fitter: " << ex.what() << std::endl;
	}

	return false;
}

#endif	// __TLIBS2_USE_MINUIT__
// ----------------------------------------------------------------------------


}

#endif
