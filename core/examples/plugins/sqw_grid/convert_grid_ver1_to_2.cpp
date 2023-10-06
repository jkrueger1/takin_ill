/**
 * converts version 1 grid data to version 2
 * @author Tobias Weber <tweber@ill.fr>
 * @date 05-jan-20
 * @license GPLv2
 */

/*
 * Grid version 2 format
 *
 * Header format:
 *     8 bytes (std::size_t): offset of index block
 *     3*8 bytes (double): data dimensions: hmin, hmax, hstep
 *     3*8 bytes (double): data dimensions: kmin, kmax, kstep
 *     3*8 bytes (double): data dimensions: lmin, lmax, lstep
 *     x bytes: metadata header
 *
 * <data block>
 * <index block>
 *
 * Data block format:
 *     repeat for each wave vector Q:
 *         4 bytes (unsigned int): number of dispersion branches
 *             repeat (number of branches times):
 *                 8 bytes (double): energy
 *                 8 bytes (double): dynamical structure factor
 *
 * Index block format:
 *     repeat for each (h,k,l) coordinate:
 *         8 bytes (std::size_t): offset into data block
 */

#define HAS_POLARISATION_DATA 0


#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>


using t_float_src = double;
using t_float_dst = double;
using t_idx = std::size_t;         // file offsets
using t_branchidx = unsigned int;  // dispersion branch indices


int main()
{
	const t_float_dst eps = 1e-8;


	std::string filenameIdx = "grid_ver1.idx";
	std::string filenameDat = "grid_ver1.bin";

	std::string filenameNewDat = "grid_ver2.bin";


	// dimensions of version 1 grid (TODO: load from config file)
	t_float_dst hmin = -0.096;
	t_float_dst hmax = 0.096;
	t_float_dst hstep = 0.002;

	t_float_dst kmin = -0.096;
	t_float_dst kmax = 0.096;
	t_float_dst kstep = 0.002;

	t_float_dst lmin = -0.096;
	t_float_dst lmax = 0.096;
	t_float_dst lstep = 0.002;


	std::cout << "Converting data file ..." << std::endl;

	std::ifstream ifDat(filenameDat);
	std::ofstream ofDat(filenameNewDat);

	std::unordered_map<t_idx, t_idx> mapIndices;

	std::size_t removedBranches = 0;


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
	ofDat << "takin_grid_data_ver2|title:TEST|author:tweber@ill.fr|date:5/feb/2020";


	while(true)
	{
		t_idx idx = ifDat.tellg();
		t_idx idxnew = ofDat.tellp();

		t_branchidx numBranches = 0;
		ifDat.read((char*)&numBranches, sizeof(numBranches));
		if(ifDat.gcount() != sizeof(numBranches) || ifDat.eof())
			break;

		mapIndices.insert(std::make_pair(idx, idxnew));
		//std::cout << "index " << std::hex << idx << " -> " << idxnew << "\n";


		// write placeholder
		t_branchidx numNewBranches = 0;
		ofDat.write((char*)&numNewBranches, sizeof(numNewBranches));


		for(t_branchidx branch=0; branch<numBranches; ++branch)
		{
#if HAS_POLARISATION_DATA != 0
			t_float_src vals[4] = { 0, 0, 0, 0 };
#else
			t_float_src vals[2] = { 0, 0 };
#endif
			ifDat.read((char*)vals, sizeof(vals));

			t_float_dst E = vals[0];
			if(std::abs(E) < eps)
				E = t_float_dst(0);

#if HAS_POLARISATION_DATA != 0
			t_float_dst w = std::abs(vals[1]) + std::abs(vals[2]) + std::abs(vals[3]);
#else
			t_float_dst w = std::abs(vals[1]);
#endif

			if(w >= eps)
			{
				t_float_dst newvals[2] = { E, w };
				ofDat.write((char*)newvals, sizeof(newvals));

				++numNewBranches;
			}
			else
			{
				++removedBranches;
			}
		}


		// seek back and write real number of branches
		t_idx lastIdx = ofDat.tellp();
		ofDat.seekp(idxnew, std::ios_base::beg);
		ofDat.write((char*)&numNewBranches, sizeof(numNewBranches));
		ofDat.seekp(lastIdx, std::ios_base::beg);
	}

	ifDat.close();

	// update index at beginning
	idx_offs = ofDat.tellp();
	ofDat.seekp(0, std::ios_base::beg);
	ofDat.write((char*)&idx_offs, sizeof(idx_offs));
	ofDat.seekp(idx_offs, std::ios_base::beg);


	std::cout << removedBranches << " branches removed (weight < eps)." << std::endl;


	std::cout << "\nConverting index file ..." << std::endl;

	std::ifstream ifIdx(filenameIdx);
	// write index block in the same file as data
	std::ofstream &ofIdx = ofDat;


	while(true)
	{
		t_idx idx = 0;

		ifIdx.read((char*)&idx, sizeof(idx));
		if(ifIdx.gcount() != sizeof(idx) || ifIdx.eof())
			break;

		auto iter = mapIndices.find(idx);
		if(iter == mapIndices.end())
		{
			std::cerr << "Error: Index " << std::hex << idx << " was not found." << std::endl;
			continue;
		}

		t_idx newidx = iter->second;
		ofIdx.write((char*)&newidx, sizeof(newidx));
	}


	return 0;
}
