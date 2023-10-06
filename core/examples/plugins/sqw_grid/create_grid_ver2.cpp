/**
 * creates a grid version 2 data file from a text data file
 * @author Tobias Weber <tweber@ill.fr>
 * @date 05-jan-20
 * @license GPLv2
 *
 * g++ -std=c++14 -I../../.. -o create_grid_ver2 create_grid_ver2.cpp
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "tlibs/string/string.h"


using t_float = double;
using t_idx = std::size_t;         // file offsets
using t_branchidx = unsigned int;  // dispersion branch indices


int main(int argc, char** argv)
{
	const t_float eps = 1e-8;

	if(argc < 12)
	{
		std::cout << "\nUsage: " << argv[0]
			<< " <input file> <output file>"
			<< " <hmin> <hmax> <hstep>"
			<< "<kmin> <kmax> <kstep>"
			<< " <lmin> <lmax> <lstep>\n"
			<< std::endl;
		return -1;
	}


	// input and output files
	std::string filename_in = argv[1];
	std::string filename_out = argv[2];

	std::ifstream ifDat(filename_in);
	std::ofstream ofDat(filename_out);


	// grid dimensions
	t_float hmin = tl::str_to_var<t_float>(std::string(argv[3]));
	t_float hmax = tl::str_to_var<t_float>(std::string(argv[4]));
	t_float hstep = tl::str_to_var<t_float>(std::string(argv[5]));

	t_float kmin = tl::str_to_var<t_float>(std::string(argv[6]));
	t_float kmax = tl::str_to_var<t_float>(std::string(argv[7]));
	t_float kstep = tl::str_to_var<t_float>(std::string(argv[8]));

	t_float lmin = tl::str_to_var<t_float>(std::string(argv[9]));
	t_float lmax = tl::str_to_var<t_float>(std::string(argv[10]));
	t_float lstep = tl::str_to_var<t_float>(std::string(argv[11]));


	// --------------------------------------------------------------------
	std::cout << "Writing header block ..." << std::endl;

	// write a dummy index file offset at the beginning (to be filled in later)
	t_idx idx_offs = 0;
	ofDat.write((char*)&idx_offs, sizeof(idx_offs));

	// write data dimensions
	ofDat.write((char*)&hmin, sizeof(hmin));
	ofDat.write((char*)&hmax, sizeof(hmax));
	ofDat.write((char*)&hstep, sizeof(hstep));

	ofDat.write((char*)&kmin, sizeof(kmin));
	ofDat.write((char*)&kmax, sizeof(kmax));
	ofDat.write((char*)&kstep, sizeof(kstep));

	ofDat.write((char*)&lmin, sizeof(lmin));
	ofDat.write((char*)&lmax, sizeof(lmax));
	ofDat.write((char*)&lstep, sizeof(lstep));

	// header metadata
	ofDat << "takin_grid_data_ver2";

	std::cout << "Grid h extents: " << hmin << " .. " << hmax << ", stepping: " << hstep << "." << std::endl;
	std::cout << "Grid k extents: " << kmin << " .. " << kmax << ", stepping: " << kstep << "." << std::endl;
	std::cout << "Grid l extents: " << lmin << " .. " << lmax << ", stepping: " << lstep << "." << std::endl;
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	std::cout << "\nWriting dispersion data block ..." << std::endl;

	std::vector<t_idx> indices;
	std::size_t removedBranches = 0;
	std::size_t totalBranches = 0;
	std::string line;

	while(std::getline(ifDat, line))
	{
		std::vector<t_float> vecHKLES;
		tl::get_tokens<t_float>(line, std::string(" \t"), vecHKLES);

		// index to this Q entry
		t_idx idxnew = ofDat.tellp();
		indices.push_back(idxnew);

		// input format: h k l E1 w1 E2 w2 ...
		t_branchidx numBranches = t_branchidx(vecHKLES.size() - 3) / 2;

		// write placeholder
		t_branchidx numNewBranches = 0;
		ofDat.write((char*)&numNewBranches, sizeof(numNewBranches));

		// TODO: check (h, k, l)
		/*t_float h = vecHKLES[0];
		t_float k = vecHKLES[1];
		t_float l = vecHKLES[2];*/

		for(t_branchidx branch=0; branch<numBranches; ++branch)
		{
			t_float E = vecHKLES[3 + branch*2 + 0];
			if(std::abs(E) < eps)
				E = t_float(0);

			t_float w = std::abs(vecHKLES[3 + branch*2 + 1]);

			if(w >= eps)
			{
				ofDat.write((char*)&E, sizeof(t_float));
				ofDat.write((char*)&w, sizeof(t_float));

				++numNewBranches;
			}
			else
			{
				++removedBranches;
			}
		}

		totalBranches += numNewBranches;

		// seek back and write real number of branches
		t_idx lastIdx = ofDat.tellp();
		ofDat.seekp(idxnew, std::ios_base::beg);
		ofDat.write((char*)&numNewBranches, sizeof(numNewBranches));
		ofDat.seekp(lastIdx, std::ios_base::beg);
	}

	ifDat.close();
	std::cout << totalBranches << " dispersion branches written." << std::endl;
	std::cout << removedBranches << " dispersion branches removed (weight < eps)." << std::endl;
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	// update index at beginning
	idx_offs = ofDat.tellp();
	ofDat.seekp(0, std::ios_base::beg);
	ofDat.write((char*)&idx_offs, sizeof(idx_offs));
	ofDat.seekp(idx_offs, std::ios_base::beg);
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	// write Q indices
	std::cout << "\nWriting Q index block ..." << std::endl;
	ofDat.write((char*)indices.data(), sizeof(t_idx)*indices.size());
	std::cout << indices.size() << " Q positions written." << std::endl;
	// --------------------------------------------------------------------

	return 0;
}
