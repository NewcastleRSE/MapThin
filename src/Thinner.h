/************************************************************************
 * MapThin version 1.11
 * Copyright 2011-2014,
 * Richard Howey
 * Institute of Genetic Medicine, Newcastle University
 *
 * richard.howey@ncl.ac.uk
 * http://www.staff.ncl.ac.uk/richard.howey/
 *
 * This file is part of MapThin, a program to thin map files.
 *
 * MapThin is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MapThin is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MapThin.  If not, see <http://www.gnu.org/licenses/>.
 ************************************************************************/


#ifndef __MAPTHIN
#define __MAPTHIN

#include <list>
#include <string>
#include <map>
#include <iostream>
#include <ostream>
#include <fstream>

using namespace std;

//! Class to store data for one SNP
struct SNP
{
	double geneticDistance; // in cM
	bool include; //include in thinned SNP file

	SNP(double & gd) : geneticDistance(gd), include(false) {};

	~SNP() {};
};

//! Class for thinning a map file.
class MapThinner
{
private:
	string filename;
	string outputFileName;
	double snpsPerCM;
	bool useBasePairPosition;
	unsigned int noMissing;
	unsigned int totalNoSNPs;
	double totalCM;
	bool writeThinnedFile;
	bool bim;
	bool search;
	bool foundUnorderedSNP;
	bool nameOnly;

	list<SNP *> theSNPs;
	map<string, list<double> > geneticDistancesInFinal; //chromosome, genetic distances, used to calc stats of final file

	ifstream readMap;
	ofstream writeMap;
	
	ofstream writeMissing; //SNPs with missing genetic distance (or base pair position)
	
public:

	MapThinner(string & fn, string & ofn, double & spc, bool & ubp, bool & no) :
	  filename(fn), outputFileName(ofn), snpsPerCM(spc), useBasePairPosition(ubp), noMissing(0), totalNoSNPs(0), totalCM(0), writeThinnedFile(true), bim(false), search(false), foundUnorderedSNP(false), nameOnly(no)
	  {
		    setBim();
			setTotalNoSNPs();
	  };

	
	~MapThinner()
	{
		for(list<SNP *>::iterator i = theSNPs.begin(); i != theSNPs.end(); ++i)
		{
			delete *i;
		};
	};

	void thin();
	void thinSNPs(string & chromosomeToThin);
	void includeSNPsForFinal();
	void outputMissing(string & chromosome, string & snpIdentifier, string & geneticDistance, string & basePairPosition, string & alleleName1, string & alleleName2);
	void displayFinalFileStats();
	void displayMissingDataStats();
	void displayWarningUnordered();
	void setSNPsPerCMFromTotalSNPs(unsigned int & totalSNPsToKeep);
	void setTotalNoSNPs();
	void setBim();
	void readLineData(ifstream & readMapFile, string & chromosome, string & snpIdentifier, string & geneticDistance, string & basePairPosition, string & alleleName1, string & alleleName2);
	void writeLineData(ofstream & writeMapFile, string & chromosome, string & snpIdentifier, string & geneticDistance, string & basePairPosition, string & alleleName1, string & alleleName2);
	unsigned int getTotalNoThinnedSNPs();
	void thinToTargetNoSNPs(unsigned int & targetThinnedSNPs);
	void thinToTargetPercentNoSNPs(double & percentToKeep);
	unsigned int tryThinning(double & snpsPerCMToTry);
	pair<double, double> getSNPsPerCMInterval(unsigned int & targetThinnedSNPs);
};

#endif

