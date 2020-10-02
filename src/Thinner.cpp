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

#include "Thinner.h"
#include "main.h"

#include <set>
#include <string>
#include <iostream>
#include <ostream>
#include <fstream>
#include <math.h>

using namespace std;

//! Thin the SNPs.
void MapThinner::thin()
{
	string chromosome, snpIdentifier, geneticDistance, basePairPosition;
	string alleleName1, alleleName2;

	double geneDis;
	double prevGeneDis = -1;
	noMissing = 0;

	readMap.open(filename.c_str());
		
	if(!readMap.is_open())
	{
		cerr<<"Cannot read map file: "<<filename << "!\n";
		exit(1);
	};

	//ensure lists used during thinning are empty
	for(list<SNP *>::iterator i = theSNPs.begin(); i != theSNPs.end(); ++i)	delete *i;
	theSNPs.clear();
	geneticDistancesInFinal.clear();

	//get first chromosome
	string prevChromosome;
	ifstream readMapFirstChromo;
	readMapFirstChromo.open(filename.c_str());
	readMapFirstChromo >> prevChromosome;
	readMapFirstChromo.close();

	if(writeThinnedFile) writeMap.open(outputFileName.c_str());

	SNP * aSNP;

	//read in map data
	do{

		readLineData(readMap, chromosome, snpIdentifier, geneticDistance, basePairPosition, alleleName1, alleleName2);

		if(readMap.eof()) break;

		//process the previous chromosome and then move onto the next
		if(chromosome != prevChromosome)
		{
			thinSNPs(prevChromosome);
			prevGeneDis = -1;
		};

		if(useBasePairPosition) geneDis = atof(basePairPosition.c_str());
		else geneDis = atof(geneticDistance.c_str());

		aSNP = new SNP(geneDis);
		if(geneDis == 0)
		{
			aSNP->include = false;
			if(writeThinnedFile) outputMissing(chromosome, snpIdentifier, geneticDistance, basePairPosition, alleleName1, alleleName2);					
			noMissing++;
		}
		else
		{
			if(geneDis < prevGeneDis && chromosome == prevChromosome)
			{
				foundUnorderedSNP = true;					
			};
			prevGeneDis = geneDis;
		};

		theSNPs.push_back(aSNP);

		prevChromosome = chromosome;

	}while(!readMap.eof());

	thinSNPs(chromosome);

	if(writeThinnedFile) writeMap.close();
	readMap.close();
	if(writeThinnedFile && noMissing > 0) writeMissing.close();

	if(outputToScreen && writeThinnedFile) displayFinalFileStats();
};

//! Thin SNPs and remove from the list.
void MapThinner::thinSNPs(string & chromosomeToThin)
{
	string chromosome, snpIdentifier, geneticDistance, basePairPosition;
	string alleleName1, alleleName2;

	ifstream readMap2;
	readMap2.open(filename.c_str());

	if(!readMap2.is_open())
	{
		cerr<<"Cannot read map file: "<<filename << "!\n";
		exit(1);
	};

	//set which SNPs are to be included in the final file
	includeSNPsForFinal();

	//denotes the start of a new chromosome for calculating the final stats
	list<double> aListGeneDis;

	//loop thro' map file until at the start of the chromosome to thin current batch of SNPs
	do{

		readLineData(readMap2, chromosome, snpIdentifier, geneticDistance, basePairPosition, alleleName1, alleleName2);

		if(chromosome == chromosomeToThin) break;

	}while(!readMap2.eof());

	//thin SNPs
	for(list<SNP *>::iterator i = theSNPs.begin(); i != theSNPs.end(); ++i)
	{
		//write SNP if it is marked to be included in the final file
		if((*i)->include)
		{
			 if(writeThinnedFile)
			 {					
					writeLineData(writeMap, chromosome, snpIdentifier, geneticDistance, basePairPosition, alleleName1, alleleName2);
			 };

			 aListGeneDis.push_back((*i)->geneticDistance);
		};

		readLineData(readMap2, chromosome, snpIdentifier, geneticDistance, basePairPosition, alleleName1, alleleName2);

		if(readMap2.eof()) break;
	};

	//add distances for this chromosome to list of genetic distances
	geneticDistancesInFinal[chromosomeToThin] = aListGeneDis;

	//delete and then empty the list
	for(list<SNP *>::iterator i = theSNPs.begin(); i != theSNPs.end(); ++i)
	{
		delete *i;
	};

	theSNPs.clear();
};

//! Outputs SNPs with missing genetic distance.
void MapThinner::outputMissing(string & chromosome, string & snpIdentifier, string & geneticDistance, string & basePairPosition, string & alleleName1, string & alleleName2)
{
	if(noMissing == 0)
	{
		string missingFileName = "missingGeneticDis.txt";
		if(useBasePairPosition) missingFileName = "missingBasePairPosition.txt";
		writeMissing.open(missingFileName.c_str());
	};

	writeLineData(writeMissing, chromosome, snpIdentifier, geneticDistance, basePairPosition, alleleName1, alleleName2);
};

//! Marks which SNPs are to be included in final file 
void MapThinner::includeSNPsForFinal()
{
	//cout << snpsPerCM << " p final\n";
	double geneDisStep = 1.0/snpsPerCM;

	if(useBasePairPosition) geneDisStep *= 1000000; //base pair position step is per 1000000 

	double marker = (*theSNPs.begin())->geneticDistance + geneDisStep;
	double prevIncludeGeneDis = (*theSNPs.begin())->geneticDistance;

	list<SNP *>::iterator i = theSNPs.begin();
	
	//include the first SNP
	(*i)->include = true;
	 ++i;

	SNP * prevSNP = *i;

	do{

		//pick a SNP to include
		if((*i)->geneticDistance > marker)
		{
			//pick closest SNP to marker or second if first is already chosen
			if(((*i)->geneticDistance - marker) <  (marker - prevSNP->geneticDistance) || prevSNP->include )
			{
				if((*i)->geneticDistance != prevIncludeGeneDis) (*i)->include = true;
				prevIncludeGeneDis = (*i)->geneticDistance;
			}
			else
			{
				if(prevSNP->geneticDistance != prevIncludeGeneDis) prevSNP->include = true;
				prevIncludeGeneDis = prevSNP->geneticDistance;
			};

			//move on marker past the geneDis of the last included SNP
			do{ marker += geneDisStep; }while(marker <= prevIncludeGeneDis);
		};

		prevSNP = *i;
		++i;
	}while(i != theSNPs.end());


};

//! Returns the final number of thinned SNPs after a SNP thinning has been done.
unsigned int MapThinner::getTotalNoThinnedSNPs()
{
	unsigned int noSNPs = 0;
	for(map<string, list<double> >::const_iterator c = geneticDistancesInFinal.begin();  c != geneticDistancesInFinal.end(); ++c)
	{
		for(list<double>::const_iterator g = c->second.begin(); g != c->second.end(); ++g)
		{
			noSNPs++;
		};
	};

	return noSNPs;
};

//! Display stats SNPs with missing genetic distances or base pair positions.
void MapThinner::displayMissingDataStats()
{

	if(noMissing > 0)
	{
			if(useBasePairPosition)
			{
				cout << "Number of SNPs with missing base pair positions: "<<noMissing<<"\n";
				cout <<"\t(Written to file missingBasePairPosition.txt)\n\n";
			}
			else
			{
				cout << "Number of SNPs with missing genetic distances: "<<noMissing<<"\n";
				cout <<"\t(Written to file missingGeneticDis.txt)\n\n";
			};
	};

};

//! Display warning of unordered SNPs on genetic distances or base pair positions.
void MapThinner::displayWarningUnordered()
{
	if(foundUnorderedSNP)
	{
		if(useBasePairPosition) cout << "Warning: SNPs in \""<<filename<<"\" are not ordered on the base pair position!\n\n";
		else cout << "Warning: SNPs in \""<<filename<<"\" are not ordered on the genetic distance!\n\n";
	};
};

//! Display stats of the SNPs in the final file.
void MapThinner::displayFinalFileStats()
{
	//get total number of SNPs
	unsigned int noSNPs = getTotalNoThinnedSNPs();
	
	cout << "Statistics: \n"
		 << "Total number of SNPs in original file: "<<totalNoSNPs<<"\n"
	     << "Number of SNPs in thinned file: "<<noSNPs<<" ("<< (100*((double)noSNPs)/((double)totalNoSNPs)) <<"%)\n";

	displayMissingDataStats();

	if(noSNPs < 2 )
	{
			cout << "\nThat's a bit too thin!\n";
			return;
	};

	double minDis = 0;
	double maxDis = 0;
	double mean = 0;
	double stdev = 0;
	double prevGD, diff;
	
	//set upper bound for minimum difference
	for(map<string, list<double> >::const_iterator c = geneticDistancesInFinal.begin();  c != geneticDistancesInFinal.end(); ++c)
	{
		minDis += *(c->second.rbegin());
	};

	//now calculate the mean
	for(map<string, list<double> >::const_iterator ch = geneticDistancesInFinal.begin();  ch != geneticDistancesInFinal.end(); ++ch)
	{
		//add stats from each chromosome in turn
		list<double>::const_iterator gd = ch->second.begin();

		prevGD = *(ch->second.begin());
		++gd;

		while(gd != ch->second.end())
		{
			diff = *gd - prevGD;
			mean += diff;

			if(diff < minDis) minDis = diff;
			if(diff > maxDis) maxDis = diff;

			prevGD = *gd;
			++gd;
		};
	};

	mean = mean/((double)(noSNPs - 1));

	//now calculate the standard dev
	for(map<string, list<double> >::const_iterator ch2 = geneticDistancesInFinal.begin();  ch2 != geneticDistancesInFinal.end(); ++ch2)
	{
		//add stats from each chromosome in turn
		list<double>::const_iterator gd2 = ch2->second.begin();

		prevGD = *(ch2->second.begin());
		++gd2;

		while(gd2 != ch2->second.end())
		{
			diff = *gd2 - prevGD;
			stdev += (diff - mean)*(diff - mean);
		
			prevGD = *gd2;
			++gd2;
		};
	};

	stdev = sqrt(stdev/((double)(noSNPs - 1)));

	cout << "\n";

	if(useBasePairPosition)
	{
		if(search) cout << "SNPs per 10^6 base pair position (in file): " << snpsPerCM << "\n";
		cout << "Mean base pair position (in file) between SNPs: "<<mean<<" bpp\n"
			 << "St. dev. of base pair position (in file) between SNPs: "<<stdev<<" bpp\n"
			 << "Range of base pair position (in file) between SNPs: ("<<minDis<<", "<<maxDis<<")\n\n";
	}
	else
	{
		if(search) cout << "SNPs per cM: " << snpsPerCM << "\n";
		cout << "Mean genetic distance between SNPs: "<<mean<<" cM\n"
			 << "St. dev. of genetic distance between SNPs: "<<stdev<<" cM\n"
			 << "Range of genetic distances between SNPs: ("<<minDis<<", "<<maxDis<<")\n\n";
	};

	displayWarningUnordered();
	
};

//! Read a line of SNP map file data
void MapThinner::readLineData(ifstream & readMapFile, string & chromosome, string & snpIdentifier, string & geneticDistance, string & basePairPosition, string & alleleName1, string & alleleName2)
{
	if(bim)
	{
		readMapFile >> chromosome >> snpIdentifier >> geneticDistance >> basePairPosition >> alleleName1 >> alleleName2;
	}
	else
	{
		readMapFile >> chromosome >> snpIdentifier >> geneticDistance >> basePairPosition;			
	};
};

//! Write a line of SNP map file data
void MapThinner::writeLineData(ofstream & writeMapFile, string & chromosome, string & snpIdentifier, string & geneticDistance, string & basePairPosition, string & alleleName1, string & alleleName2)
{
	if(nameOnly)
	{
		writeMapFile << snpIdentifier <<"\n";	
	}
	else
	{
		if(bim)
		{
			writeMapFile  << chromosome << "\t" << snpIdentifier << "\t" << geneticDistance << "\t" << basePairPosition <<"\t" << alleleName1 <<"\t"<< alleleName2<<"\n";;
		}
		else
		{
			writeMapFile  << chromosome << "\t" << snpIdentifier << "\t" << geneticDistance << "\t" << basePairPosition <<"\n";		
		};	
	};
};

//! Sets whether the map file is a .bim file or not (used with binary files)
void MapThinner::setBim()
{
	unsigned int length = filename.length();
	string fileExtension;

	if(length >= 4) fileExtension = filename.substr(length-4,4);
	else return;

	//determine the type of pedigree file
	if(fileExtension[0] == '.' &&
		(fileExtension[1] == 'b' || fileExtension[1] == 'B') &&
		(fileExtension[2] == 'i' || fileExtension[2] == 'I') &&
		(fileExtension[3] == 'm' || fileExtension[3] == 'M'))
	{
		bim = true;
	};

};

//! Sets the total number of SNPs in original map file and the total cM distance
void MapThinner::setTotalNoSNPs()
{
	string chromosome, snpIdentifier, geneticDistance, basePairPosition;
	string alleleName1, alleleName2;
	
	ifstream readMap3;
	readMap3.open(filename.c_str());

	if(!readMap3.is_open())
	{
		cerr<<"Cannot read map file: "<<filename << "!\n";
		exit(1);
	};

	//get first chromosome
	string prevChromosome, prevGeneticDistance;
	ifstream readMapFirstChromo;
	readMapFirstChromo.open(filename.c_str());
	readMapFirstChromo >> prevChromosome;
	readMapFirstChromo.close();

	do{

		readLineData(readMap3, chromosome, snpIdentifier, geneticDistance, basePairPosition, alleleName1, alleleName2);

		if(readMap3.eof()) break;

		if(chromosome != prevChromosome)
		{
			if(useBasePairPosition) totalCM += atof(basePairPosition.c_str());
			else totalCM += atof(prevGeneticDistance.c_str());
		};

		totalNoSNPs++;

		prevChromosome = chromosome;
		prevGeneticDistance = geneticDistance;

	}while(!readMap3.eof());

};

//! Sets SNPs per cM based on the total no. of SNPs to keep
void MapThinner::setSNPsPerCMFromTotalSNPs(unsigned int & totalSNPsToKeep)
{
	if(totalCM != 0) snpsPerCM = (double)(totalSNPsToKeep)/totalCM;
};


//! Use a bisection search to thin SNPs to the total required
void MapThinner::thinToTargetNoSNPs(unsigned int & targetThinnedSNPs)
{
	search = true;
	if(!(targetThinnedSNPs < totalNoSNPs && targetThinnedSNPs > 0))
	{
		cerr << "The number of SNPs to keep must be between 0 and "<<totalNoSNPs<<"!\n";
		exit(1);
	};

	setSNPsPerCMFromTotalSNPs(targetThinnedSNPs);
	pair<double, double> snpsPerCMInterval = getSNPsPerCMInterval(targetThinnedSNPs);

	if((snpsPerCMInterval.first == -1 || snpsPerCMInterval.first == 0) && snpsPerCMInterval.second == 0)
	{
		cout << "Statistics: \n"
		     << "Total number of SNPs in original file: "<<totalNoSNPs<<"\n\n";
		displayMissingDataStats();
		displayWarningUnordered();

		cerr << "Failed to thin SNPs for these settings!\n\n";
		if(snpsPerCMInterval.first == -1)
		{
			if(useBasePairPosition) cerr << "Consider using genetic distance instead (do not use the -b option)!\n\n";
			else cerr << "Consider using the base pair position option (-b)!\n\n";
		};
		exit(1);
	};

	double trySNPsPerCM;
	unsigned int nextGuess;
	unsigned int count = 0;

	do{		
		//try mid way point as next guess
		trySNPsPerCM = (snpsPerCMInterval.first + snpsPerCMInterval.second)*0.5;

		nextGuess = tryThinning(trySNPsPerCM);

		if(nextGuess == targetThinnedSNPs || (snpsPerCMInterval.second - snpsPerCMInterval.first) < 1e-6) break;
		else if(nextGuess > targetThinnedSNPs) snpsPerCMInterval.second = trySNPsPerCM;
		else snpsPerCMInterval.first = trySNPsPerCM;

		count++;

	}while(count < 100);

	//thin SNPs and write to file using the best found SNPs per cM to achieve target no of SNPs
	writeThinnedFile = true;
	thin();
};

//! Thins to a target percentage of SNPs
void MapThinner::thinToTargetPercentNoSNPs(double & percentToKeep)
{
	unsigned int targetNoSNPs = (unsigned int)((double)(totalNoSNPs)*(percentToKeep*0.01) + 0.5);
	thinToTargetNoSNPs(targetNoSNPs);
};

//! Get an interval for the target number of SNPs
pair<double, double> MapThinner::getSNPsPerCMInterval(unsigned int & targetThinnedSNPs)
{
	unsigned int initialGuess = tryThinning(snpsPerCM);

	if(targetThinnedSNPs > totalNoSNPs - noMissing) return make_pair(-1, 0);

	unsigned int otherGuess;
	double aBoundsnpsPerCM = snpsPerCM;
	double otherBound;
	unsigned int count = 0;
	double step = 0.5;
	if(initialGuess > targetThinnedSNPs) step = -0.5;

	otherBound = aBoundsnpsPerCM;
	
	do{
		otherBound += step;
		if(otherBound < 0)
		{
			otherBound = (otherBound - step)*0.5;
		};

		otherGuess = tryThinning(otherBound);

		count++;

		if(count > 100) return make_pair(0, 0);//failed to find a bounding interval

	}while(!((otherGuess < targetThinnedSNPs && initialGuess > targetThinnedSNPs)
		      || (otherGuess > targetThinnedSNPs && initialGuess <= targetThinnedSNPs)));

	if(otherGuess <= targetThinnedSNPs && initialGuess > targetThinnedSNPs) return make_pair(otherBound, aBoundsnpsPerCM);
	else return make_pair(aBoundsnpsPerCM, otherBound);
};

//! Trys to thin the SNPs and returns the total number of SNPs in the thinned SNP file.
unsigned int MapThinner::tryThinning(double & snpsPerCMToTry)
{
	//set variables to try a thinning
	writeThinnedFile = false;

	//try a thinning with with value for SNP per cM
	snpsPerCM = snpsPerCMToTry;
	thin();

	return getTotalNoThinnedSNPs();
};


