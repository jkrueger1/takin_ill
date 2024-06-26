/**
 * mini scan monitor
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 2015
 * @license GPLv2
 *
 * ----------------------------------------------------------------------------
 * Takin (inelastic neutron scattering software package)
 * Copyright (C) 2017-2023  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2013-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; version 2 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 * ----------------------------------------------------------------------------
 */

#include <iostream>
#include <iomanip>
#include <fstream>

#include "../../tlibs/net/tcp.h"
#include "../../tlibs/log/log.h"
#include "../../tlibs/string/string.h"
#include "dialog.h"


using namespace tl;
typedef float t_real;

static t_real dCtr = 0.;
static t_real dMon = 0.;   // or time
static t_real dSel = 0.;
static bool bCountToMon = 1;


static const std::string strCtr = "nicos/ctr1/value";
static const std::string strTim = "nicos/timer/value";
static const std::string strSel = "nicos/timer/preselection";


void refresh()
{
	t_real dProgress = 0.;
	t_real dExpCtr = 0.;

	if(!float_equal(dSel, t_real(0.)))
	{
		dProgress = dMon / dSel;
		if(!float_equal(dProgress, t_real(0.)))
			dExpCtr = dCtr / dProgress;
	}

	if(dProgress < t_real(0.)) dProgress = t_real(0.);
	else if(dProgress > t_real(1.)) dProgress = t_real(1.);


	std::ostringstream ostr;
	ostr.precision(2);

	std::string strMon = "Monitor:  ";
	if(!bCountToMon)
		strMon = "Time:     ";
	ostr << "Counts:   " << std::fixed << dCtr << " +- " << std::sqrt(dCtr) << "\n";
	ostr << "Expected: " << std::fixed << int(std::round(dExpCtr)) << " (" << dExpCtr << ")\n";
	ostr << strMon << std::fixed << dMon << " of " << dSel
		<< (bCountToMon ? " counts" : " seconds") << "\n";
	ostr << "Progress: " << std::fixed << dProgress*t_real(100.) << " %";
	set_progress(int(std::round(dProgress*t_real(100.))), ostr.str());
}


void disconnected(const std::string& strHost, const std::string& strSrv)
{
	//log_info("Disconnected from ", strHost, " on port ", strSrv, ".");;
}


void connected(const std::string& strHost, const std::string& strSrv)
{
	//log_info("Connected to ", strHost, " on port ", strSrv, ".");
}


void received_sics(const std::string& strMsg)
{
	std::pair<std::string, std::string> pair = split_first<std::string>(strMsg, "=", true);

	if(str_is_equal<std::string>(pair.first, "counter.Monitor 1"))
		dMon = str_to_var<t_real>(pair.second);
	else if(str_is_equal<std::string>(pair.first, "counter.Counts"))
		dCtr = str_to_var<t_real>(pair.second);
	else if(str_is_equal<std::string>(pair.first, "counter.Preset"))
		dSel = str_to_var<t_real>(pair.second);

	refresh();
}


void received_nicos(const std::string& strMsg)
{
	std::pair<std::string, std::string> pair = split_first<std::string>(strMsg, "=", true);
	std::string& strVal = pair.second;
	std::size_t iBeg = strVal.find_first_of("0123456789e.+-");
	if(iBeg == std::string::npos)
		return;
	strVal = strVal.substr(iBeg);

	if(pair.first == strTim)
		dMon = str_to_var<t_real>(strVal);
	else if(pair.first == strCtr)
		dCtr = str_to_var<t_real>(strVal);
	else if(pair.first == strSel)
		dSel = str_to_var<t_real>(strVal);

	refresh();
}


int main(int argc, char** argv)
{
	bool bUseNicos = 1;

	if(argc == 3)
	{
		bUseNicos = 1;
		bCountToMon = 0;
	}
	else if(argc == 5)
	{
		bUseNicos = 0;
		bCountToMon = 1;
	}
	else
	{
		std::cerr << "Usage: \n"
			<< "\t" << argv[0] << " <nicos server> <port>\n"
			<< "\t" << argv[0] << " <sics server> <port> <login> <password>\n\n"
			<< "\t e.g.: " << argv[0] << " mira1.mira.frm2 14869"
			<< std::endl;
		return -1;
	}


	if(!open_progress("Scan Progress"))
	{
		log_err("Cannot open progress dialog.");
		return -1;
	}


	TcpTxtClient<> client;
	if(bUseNicos)
		client.add_receiver(received_nicos);
	else
		client.add_receiver(received_sics);
	client.add_disconnect(disconnected);
	client.add_connect(connected);


	if(!client.connect(argv[1], argv[2]))
	{
		log_err("Error: Cannot connect to instrument server.");
		return -1;
	}


	if(bUseNicos)	// nicos
	{
		for(const std::string& strKey : {strSel, strTim, strCtr})
		{
			client.write(strKey + "?\n");
			client.write(strKey + ":\n");
		}

		client.wait();
	}
	else		// sics
	{
		std::string strLogin = argv[3];
		std::string strPwd = argv[4];
		client.write(strLogin + " " + strPwd + "\n");
		std::this_thread::sleep_for(std::chrono::milliseconds(250));

		while(1)
		{
			client.write("counter getcounts\ncounter getmonitor 1\ncounter getpreset\n");
			std::this_thread::sleep_for(std::chrono::milliseconds(750));
		}
	}

	close_progress();
	return 0;
}
