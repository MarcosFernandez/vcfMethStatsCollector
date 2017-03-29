/*
 * main.cpp
 *
 *  Created on: 28 Mar, 2017
 *      Author: marcos
 */
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <map>
#include <sstream>
#include <algorithm>
#include <fstream>

using namespace std;

/**
 * \brief splits a string and leaves the results in a vector
 * \param s string to be split
 * \param delim character delimitate
 * \param elems vector of string elements to be fulfilled with the different fields
 */
vector<string> split(std::string s, char delim, vector<string> elems)
{
	std::stringstream ss(s);
	std::string item;

	while(getline(ss, item, delim))
	{
		elems.push_back(item);
	}
	return elems;
}

/**
 * \brief Splits a string according to a given char separator
 * \param s string to be split
 * \param delim character delimitate
 * \returs vector of elements
 */
vector<string> split(const string &s, char delim)
{
		vector<string> elems;
	    return split(s, delim, elems);
}

/**
 * \brief Eval Status Position according to the reference and genotype
 * \param reference Reference Nucleotide
 * \param genotype Genotype Called
 * \param typeMutation, in case of snp reports the type of mutation
 * \returns 0 Equal To reference 1 Snp Mutation 2 Indel 3 Multiallelic
 */
unsigned int statusPosition(const string & reference, const string & genotype,string & typeMutation)
{
	/*1. Genotype Is Not Like reference */
	if(genotype.compare(".") != 0)
	{
		unsigned int refLen = reference.length();
		unsigned int genLen = genotype.length();
		/*2. SNP Mutation*/
		if(refLen == 1 && genLen == 1)
		{
			typeMutation = reference + ">" + genotype;
			return 1;
		}
		else if (refLen == 1 && genLen > 1)
		{
			if(genotype.find(",") != std::string::npos)
			{
				/*3. Multiallelic*/
				return 3;
			}
		}
		else if(refLen > 1 || genLen > 1)
		{
			/*4.InDel*/
			return 2;
		}
	}
	/*5. Genotype is like the reference*/
	return 0;
}


/**
 * \brief Application Usage
 */
void printHelp()
{
	cout << "vcfMethStatsCollector Collects basic variants stats from a VCF input." << endl;
	cout << "vcfMethStatsCollector -j jsonOutputfile" << endl;
	cout << "Example: bcftools file.bcf | vcfMethStatsCollector -j jsonOutputfile" << endl;
}

/**
 * \brief main function
 * \brief Reads VCF File From Standard Input and Estimates a set of stats
 */
int main(int argc, char *argv[])
{
	/*0. Variables Initialization*/
	unsigned int nSnps = 0;
	unsigned int nInDels = 0;
	unsigned int nMultiallelic = 0;
	map<unsigned int,unsigned int> coverageVariants;
	map <unsigned int,unsigned int> genotypeQualityVariants;
	map <string,unsigned int> mutationChanges;
	map <string,unsigned int> chromosomeVariants;

	string jsonFile;

	/*0. Check Arguments */
	/*0.1 Arguments checking*/
	if (argc == 1)
	{
		printHelp();
		return 1;
	}

	/*0.2 Arguments processing */
	for (int i = 1; i < argc; i++)
	{
		if (string(argv[i]).compare("-j") == 0)
		{
			jsonFile = string(argv[i + 1]);
		}
	}

	/*0.3 Checks json file */
	if(jsonFile.empty())
	{
		cout << "Sorry!! Json file was not specified!!" << endl;
		return 1;
	}

	/*1. Read From Standard Input */
	ios_base::sync_with_stdio(false); //  Don't sync C++ and C I/O

	string line;
	string typeOfMutation;
	/*1.1 Loop while there is something to read at standard input*/
	while(cin)
	{
		/*1.1.1 Get line from standard input */
	    getline(cin, line);

	    if(!line.empty())
	    {
	    	/*1.1.1.1 Process those line which are not comment*/
	    	if(line[0]!= '#')
	    	{
	    		/*1.1.1.1.1 Split input line in a vector <string> */
	    		vector<string> fields = split(line,'\t');
	    		/* Line Format O=CHROM 1=POS 2=ID 3=REF 4=ALT 5=QUAL 6=FILTER 7=INFO 8=FORMAT 9=SAMPLE*/
	    		bool updateVariantsStats = false;
	    		/*A. Process Status Position */
	    		switch (statusPosition(fields[3],fields[4],typeOfMutation))
	    		{
	    			case 1:
	    			{
	    				nSnps ++;

	    				/* A.A Mutation profiles*/
						if ( mutationChanges.find(typeOfMutation) == mutationChanges.end() )
						{
							//1. Mutation not found in the map
							mutationChanges[typeOfMutation] = 1;
						} else {
							//2. Mutation already found in the map
							mutationChanges[typeOfMutation] ++;
						}
						updateVariantsStats = true;
	    				break;
	    			}
	    			case 2:
	    			{
	    				nInDels ++;
	    				updateVariantsStats = true;
	    				break;
	    			}
	    			case 3:
	    			{
	    				nMultiallelic ++;
	    				updateVariantsStats = true;
	    				break;
	    			}
	    		}

	    		/*1.1.1.1.2 Update variants */
	    		if(updateVariantsStats)
	    		{
	    			/* A.A Compute coverage Snps*/
	    			unsigned int coverage = atoi(split(fields[9],':')[1].c_str());
	    			if(coverageVariants.find(coverage) == coverageVariants.end())
	    			{
	    				//Coverage Value not found
	    				coverageVariants[coverage] = 1;
	    			} else {
	    				//Coverage Value found
	    				coverageVariants[coverage] ++;
	    			}

	    			/* A.B Genotype Quality Snps*/
	    			unsigned int quality = atoi(fields[5].c_str());
	    			if ( genotypeQualityVariants.find(quality) == genotypeQualityVariants.end() )
	    			{
	    				//1. Mutation not found in the map
	    				genotypeQualityVariants[quality] = 1;
	    			} else {
	    				//2. Mutation already found in the map
	    				genotypeQualityVariants[quality] ++;
	    			}

	    			/* A.C  Record Chromosome changes */
	    			string chromosome = fields[0];
	    			if ( chromosomeVariants.find(chromosome) == chromosomeVariants.end() )
	    			{
	    				//1. Mutation not found in the map
	    				chromosomeVariants[chromosome] = 1;
	    			} else {
	    				//2. Mutation already found in the map
	    				chromosomeVariants[chromosome] ++;
	    			}
	    		}
	    	}
	    }
	}

	/*1.1.2 Output Results to JSON File*/
	ofstream jsonOutput(jsonFile.c_str(),std::ios_base::out);

	/*1.1.2.1 Basic Snps Stats*/
	jsonOutput << "{" << endl;

	jsonOutput << "    \"snps\":"<< nSnps << "," << endl;
	jsonOutput << "    \"indels\":"<< nInDels << "," << endl;
	jsonOutput << "    \"multialelli\":"<< nMultiallelic << "," << endl;

	/*1.1.2.2 Coverage Variants*/
	jsonOutput << "    \"coverageVariants\": {"<< endl;
	map<unsigned int,unsigned int>::iterator final_iter = coverageVariants.end();
	--final_iter;

	for (map<unsigned int,unsigned int>::iterator it=coverageVariants.begin(); it!=coverageVariants.end(); ++it)
	{
		jsonOutput << "        \"" << it->first << "\":" << it->second;
		if (it != final_iter) jsonOutput << ",";
		jsonOutput << endl;
	}

	jsonOutput << "    },"<< endl;

	/*1.1.2.3 Genotype Quality Variants*/
	jsonOutput << "    \"qualityVariants\": {"<< endl;
	final_iter = genotypeQualityVariants.end();
	--final_iter;

	for (map<unsigned int,unsigned int>::iterator it=genotypeQualityVariants.begin(); it!=genotypeQualityVariants.end(); ++it)
	{
		jsonOutput << "        \"" << it->first << "\":" << it->second;
		if (it != final_iter) jsonOutput << ",";
		jsonOutput << endl;
	}

    jsonOutput << "    },"<< endl;

	/*1.1.2.4 Mutation Changes*/
	jsonOutput << "    \"mutations\": {"<< endl;
	map <string,unsigned int>::iterator final_iterator = mutationChanges.end();
	--final_iterator;

	for (map <string,unsigned int>::iterator it=mutationChanges.begin(); it!=mutationChanges.end(); ++it)
	{
		jsonOutput << "        \"" << it->first << "\":" << it->second;
		if (it != final_iterator) jsonOutput << ",";
		jsonOutput << endl;
	}

    jsonOutput << "    },"<< endl;

    /*1.1.2.4 Chromosome Variants*/
    jsonOutput << "    \"chromosomeVariants\": {"<< endl;
    final_iterator = chromosomeVariants.end();
    --final_iterator;

    for (map <string,unsigned int>::iterator it=chromosomeVariants.begin(); it!=chromosomeVariants.end(); ++it)
    {
    	jsonOutput << "        \"" << it->first << "\":" << it->second;
    	if (it != final_iterator) jsonOutput << ",";
    	jsonOutput << endl;
    }

    jsonOutput << "    }"<< endl;


    jsonOutput << "}"<< endl;

	jsonOutput.close();


	/*1.1.3 Output Results to standard output*/

	cout << "Snps          :" << nSnps << endl;
	cout << "InDels        :" << nInDels << endl;
	cout << "Multialellic  :" << nMultiallelic << endl;


	cout << "" << endl << "Coverage" << endl;
    for (map<unsigned int,unsigned int>::iterator it=coverageVariants.begin(); it!=coverageVariants.end(); ++it)
		    cout << it->first << " => " << it->second << '\n';

    cout << "" << endl << "Genotype Quality" << endl;
	for (map <unsigned int,unsigned int>::iterator it=genotypeQualityVariants.begin(); it!=genotypeQualityVariants.end(); ++it)
		    std::cout << it->first << " => " << it->second << '\n';

	cout << "" << endl << "Mutation Changes" << endl;
	for (map <string,unsigned int>::iterator it=mutationChanges.begin(); it!=mutationChanges.end(); ++it)
		    std::cout << it->first << " => " << it->second << '\n';

	cout << "" << endl << "Chromosome" << endl;
	for (map <string,unsigned int>::iterator it=chromosomeVariants.begin(); it!=chromosomeVariants.end(); ++it)
		    std::cout << it->first << " => " << it->second << '\n';



	return 0;
}
