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


#include <iostream>
#include <ostream>
#include <set>
#include <string>


using namespace std; // initiates the "std" or "standard" namespace
 
#include "main.h"
#include "Thinner.h"

bool outputToScreen = true; 

//! Output program title to screen
void header()
{
	if(outputToScreen) cout << "\nMapThin (v1.11): Thinning your map files!\n"
		<< "-----------------------------------------------------------------\n"
		<< "Copyright 2011-2014 Richard Howey, GNU General Public License, v3\n"
		<< "Institute of Genetic Medicine, Newcastle University\n\n";
};

//! Output program usage to screen
void usage()
{
		header();
	 	
		cout << "Usage:\n\t ./mapthin [options] data-in.map data-out.map\n\n"

		<< "Options:\n"
		<< "  -t x          -- SNPs per cM, x\n"
		<< "  -s y          -- Total no. of SNPs to keep, y\n"
		<< "  -p z          -- Percentage of SNPs to keep, z\n"	
		<< "  -b [w]        -- Use base pair position [with w SNPs per 10^6 bpp in file]\n"	
		<< "  -n            -- Output the name of the SNPs only\n"	
		<< "  -so           -- suppress output to screen\n\n"
		<< "Default Options:\n"
		<< "  -t 2.4\n\n";

};

//! The start of the program
int main(int argc, char * argv[])
{
	int argcount = 1;
	string option;
	string filename = "";
	string outputFileName = "";
	double snpsPerCM = 2.4;
	unsigned int totalSNPsToKeep = 0;
	double percentToKeep = 0;
	unsigned int format = 1; // 1 = PLINK, 2 = Merlin
	bool useBasePairPosition = false; //if this is set to true then read snps per cM as base pair position thro'out
	outputToScreen = true;
	bool nameOnly = false;

	//set given options
	while(argcount < argc && argv[argcount][0] == '-')
    {
		option = argv[argcount];

		if(option ==  "-b")
		{		
			useBasePairPosition = true;
			snpsPerCM = 6.7;

			argcount++;
			option = argv[argcount];

			//check if a number for the SNPs per Mbase was specified, if not process option as before
			if(option.substr(0,1) != "-")
			{
				 if(!(option.length() >= 4 && (option.substr(option.length()-4, 4) == ".map" || option.substr(option.length()-4, 4) == ".bim" 
					 || option.substr(option.length()-4, 4) == ".MAP" || option.substr(option.length()-4, 4) == ".BIM")))
				 {
					snpsPerCM = atof(argv[argcount]);
					option = "--";
				 }
				 else break; //end of options reached

			};
		};	


		if(option ==  "-t")
		{			
			argcount++; if(argcount >= argc) break;
			snpsPerCM = atof(argv[argcount]);			
		}
		else if(option ==  "-s")
		{			
			argcount++; if(argcount >= argc) break;
			totalSNPsToKeep = atoi(argv[argcount]);			
		}
		else if(option ==  "-p")
		{			
			argcount++; if(argcount >= argc) break;
			percentToKeep = atof(argv[argcount]);			
		}		
		else if(option == "-so") outputToScreen = false;
		else if(option == "-n") nameOnly = true;
		else if(option == "--") {}
		else
		{
			header();
    		cerr << "Unrecognised command line switch: " << option << "\n\n";			
    		exit(1);
		};

		argcount++;
	};

	if(argcount < argc) filename = argv[argcount++];	
	if(argcount < argc) outputFileName = argv[argcount];	

	if(filename == "" || outputFileName == "")
	{
		usage();
		exit(0);
	};	

	
	//output options to screen
	if(outputToScreen)
	{
		header();
		cout << "Parameters:\n";
		cout << "Input file: "<<filename <<"\n";
		cout << "Output file: "<<outputFileName;
		if(nameOnly) cout << " - SNP names only\n";
		else cout << "\n";
		
		if(totalSNPsToKeep > 0) cout << "Total SNPs to keep: "<< totalSNPsToKeep <<"\n";
		else if(percentToKeep > 0) cout << "Percentage of SNPs to keep: "<< percentToKeep <<"\n";
		else if(!useBasePairPosition) cout << "SNPs per cM: "<< snpsPerCM <<"\n";
		else cout << "SNPs per 10^6 base pair position (in file): "<< snpsPerCM <<"\n";
		if(useBasePairPosition && (totalSNPsToKeep > 0 || percentToKeep > 0)) cout << "Using base pair position\n";
		cout << "\n";
	};

	//check parameters
	if(percentToKeep != 0 && !(percentToKeep < 100 && percentToKeep > 0))
	{
		cerr << "The percentage of SNPs to keep must be between 0 and 100!\n";
		exit(1);
	};

	//create mapthinner and then thin
	MapThinner mapThinner(filename, outputFileName, snpsPerCM, useBasePairPosition, nameOnly);

	if(totalSNPsToKeep > 0)
	{		
		mapThinner.thinToTargetNoSNPs(totalSNPsToKeep);
	}
	else if(percentToKeep > 0)
	{		
		mapThinner.thinToTargetPercentNoSNPs(percentToKeep);
	}
	else
		mapThinner.thin();

};

