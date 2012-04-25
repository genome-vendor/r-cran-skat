/*************************************************************
 *
 * Moving window project
 * File: setid_bim_index.h
 * Date: Dec 9, 2010
 * Author: Larissa Miropolsky
 *
 * Description:
 *   create indexing - hush table between SetID.txt file and *.bim file
 *	 [setid_row#] = bim_row#
 *
 * Ver 1.1
 **************************************************************/
#include <cstring>
#include <stdio.h>
#include <string.h>
#include <bitset>
#include <math.h>
#include <fstream>
#include <iostream>

#include "setid_bim_index.h"

//=======================================================================
// This function split "str" by "delimiters" and put the result to "tokens"
// Inputs:
// str - source string
// tokens - target vector of strings
// delimeters - string separator
// Notes:
// "tokens" - vector of strings, so to move inside use:
// tokens.at(0).c_str(); tokens.at(1).c_str().
// important clear "tokens" after each iteration: tokens.clear();
// otherwise it will append new results to existing.
//=======================================================================
void Hasht::Tokenize(const std::string& str,
                      std::vector<std::string>& tokens,
                      const std::string& delimiters )
{
    // Skip delimiters at beginning.
    std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    std::string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (std::string::npos != pos || std::string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}

//=======================================================================
// Constructor function
// Inputs:
// setID - name of file includes path(full or relative).
//		This file should contain two columns:
//		first setid, second snpid, delimited by "tab" or by "space"
// bim - name of file includes path(full or relative). "*.bim" file
// mwa - name of future file includes path(full or relative). "*.mwa" file.
//		here this "*.mwa" file not will be created and not in use
//		- it just to locate "*.LOG" file in the same place as future "*.mwa".
// myerror - allocated memory to put there the final information
//		if some error happen during the run.
//=======================================================================
Hasht::Hasht(char * setID, char * bim, char* mwa, int * myerror)
{

	*myerror = NO_ERRORS;

	this->m_setidfile = setID;

	this->m_bimfile = bim;
	this->m_num_of_snps_insetid_org = -1;

	char log_filename[1000];
	memset(log_filename,'\0',sizeof(log_filename));
	strcat(log_filename,mwa);
	strcat(log_filename,"_LOG.txt");
	this->m_log.open(log_filename);


	this->upload_snpid_from_bim(myerror);
	if (*myerror != 0)
		return;

	this->upload_snpid_from_setid_build_hash(myerror);
	if (*myerror != 0)
		return;

	for (int i = 0; i < m_num_of_snps; ++ i)
		delete [] this->m_bimf_snpsid[i];
	delete [] this->m_bimf_snpsid;
	delete [] this->m_bimf_sorted;

	this->m_log.close();

}
//=======================================================================
// Destructor function
// Free memory that was dynamically allocated during the run
//=======================================================================
Hasht::~Hasht()
{

	delete [] this->m_hash_table;
	///delete [] this->m_snp_sets;  // Don't remove this line! This delete will be applyed from bed_reader module destructor

	for (int i = 0; i < m_num_of_snps_insetid; ++ i)
		delete [] this->m_setidf_setid[i];
	delete [] this->m_setidf_setid;
}

//=========================================================================
// This function creates lookup table between "setID" file and "*.bim" file
// Reads line by line "setID" file
// For every "snpID" finds ralated line in "*.bim" file
// Inputs:
// myerror - allocated memory to put there the final information
//		if some error happen during the run.
//=========================================================================
void Hasht::upload_snpid_from_setid_build_hash(int * myerror)
{
	std::string line;
	std::vector<std::string> tokens;
	std::vector<std::string> tokens2, tokens3;

	Get_Num_of_SNPs_in_SetID(myerror);
	if(this->m_num_of_snps_insetid_org < 0){
		*myerror = CANT_OPEN_SETID_FILE4READ;
		return;
	}

	this->m_hash_table = new int[this->m_num_of_snps_insetid_org +1];
	this->m_setidf_setid = new char*[this->m_num_of_snps_insetid_org +1];

	this->m_setid.open(this->m_setidfile);
	if (!this->m_setid)
	{
		*myerror = CANT_OPEN_SETID_FILE4READ;
		return;
	}

	int index = 0;
	this->m_log << "#==== The following SNP IDs requested by SetId file not found in *.Bim file ====#" << std::endl;

	while (!this->m_setid.eof( ) )
    {
		tokens.clear();
		getline(this->m_setid, line);
		Tokenize(line, tokens, "\t");
		if (tokens.size() < 2)
		{
			tokens.clear();
			Tokenize(line, tokens, " ");
		}


		if (tokens.size() >= 2)
		{
			tokens2.clear();
			Tokenize(tokens.at(1).c_str(), tokens2, " ");
			tokens3.clear();
			Tokenize(tokens2.at(0).c_str(), tokens3, "\r");
			//if one or more spaces at the end of the line -
			//spaces at the end of the line don't be taken in account.
			//tokens2.at(0).c_str() instead tokens.at(1).c_str()
			int result_of_search = this->binsearch(tokens3.at(0).c_str());
			if (result_of_search == -1)
				this->m_log << line << std::endl;
			else
			{
				this->m_hash_table[index] = result_of_search;
				m_setidf_setid[index] = new char[SNP_ID_SIZE];
				strcpy (m_setidf_setid[index] , tokens.at(0).c_str()); //copy the setId
				index ++;
			}
		}
	}
	m_num_of_snps_insetid = index;
	this->m_setid.close();

}

//=========================================================================
// Added by SLEE
//	This function get number of SNPs in the SNP set.
//=========================================================================
void Hasht::Get_Num_of_SNPs_in_SetID(int * myerror){

	std::ifstream temp_file;
	std::string line;
	temp_file.open(this->m_setidfile);
	if (!temp_file)
	{
		*myerror =  CANT_OPEN_SETID_FILE4READ;
		return;
	}

	this->m_num_of_snps_insetid_org = 0;//0;
	while (!temp_file.eof( ) ) {
		getline(temp_file, line);
		this->m_num_of_snps_insetid_org++;
	}

	temp_file.close();
}

//=======================================================================
// This function looking for "source" - snpID from "setID" file
// inside of array of "snpID"s from "*.bim" file
// Inputs:
// source - char* , snpID from "setID" file
// Outputs:
// Index of this "source" inside of array of "snpID"s from "*.bim" file
//=======================================================================
int Hasht::binsearch(const char* source)//int ar[],int size,
{
	int lb = 0,	ub = this->m_num_of_snps-1,	mid;             //lb=>lower bound,ub=>upper bound
	for( ;lb <= ub; )
	{
		mid = (lb + ub) / 2;
		if(strcmp( this->m_bimf_snpsid[m_bimf_sorted[mid]], source) == 0)
			return m_bimf_sorted[mid];  // FOUND!!
		else if(strcmp( this->m_bimf_snpsid[m_bimf_sorted[mid]], source) > 0)
			ub = mid - 1;
		else if(strcmp( this->m_bimf_snpsid[m_bimf_sorted[mid]], source) < 0)
			lb = mid + 1;
	}
	return -1 ; //cout<<"\n SEARCH UNSUCCESSFUL";
}
//=======================================================================
// This function reads "*.bim" file,
// creates "m_bimf_snpsid" array of char* - all snpIDs from "*.bim" file
// creates "m_bimf_sorted" int array of indexes to "m_bimf_snpsid" to sort it.
// creates "m_snp_sets" array of SNP_info objects that will hold info during future "*.mwa" preparation
// of how many genotypes of every type per snp - to calculate minor and major
// Inputs:
// myerror - allocated memory to put there the final information
//		if some error happen during the run.
//=======================================================================
void Hasht::upload_snpid_from_bim(int * myerror)
{
	std::string line;
	std::vector<std::string> tokens;

	this->m_bim.open(this->m_bimfile);
	if (!this->m_bim)
	{
		*myerror = CANT_OPEN_BIM_FILE4READ;
		return;
	}

	this->m_num_of_snps = -1;//0;
	while (!this->m_bim.eof( ) )
    {
		getline(this->m_bim, line);
		this->m_num_of_snps++;
	}

	this->m_bim.close();
	//----------------------------------------------

	m_snp_sets = new SNP_info[this->m_num_of_snps];
	for (int j = 0; j < m_num_of_snps;++j)
	{
		this->m_snp_sets[j].letters[0] = NULL;
		this->m_snp_sets[j].letters[1] = NULL;
		this->m_snp_sets[j].total_counter_per_letter[0] = 0;
		this->m_snp_sets[j].total_counter_per_letter[1] = 0;
		this->m_snp_sets[j].line_counter_per_letter[0] = 0;
		this->m_snp_sets[j].line_counter_per_letter[1] = 0;

	}
	//----------------------------------------------

	this->m_bim.open(this->m_bimfile);
	this->m_bim.seekg (0, std::ios::beg);
	if (!this->m_bim)
	{
		*myerror = CANT_OPEN_BIM_FILE4READ;
		return;
	}


	m_bimf_snpsid = new char*[this->m_num_of_snps];
	m_bimf_sorted = new int[this->m_num_of_snps];

	for (int i = 0; i < m_num_of_snps;++i)
    {
		tokens.clear();
		getline(this->m_bim, line);
		Tokenize(line, tokens, "	");
		strcpy (this->m_snp_sets[i].snp_id , tokens.at(1).c_str());
		this->m_snp_sets[i].letters[0] = tokens.at(4).c_str()[0];
		this->m_snp_sets[i].letters[1] = tokens.at(5).c_str()[0];
		this->m_snp_sets[i].flag = 0;

		m_bimf_snpsid[i] = new char[SNP_ID_SIZE];
		strcpy (m_bimf_snpsid[i] , tokens.at(1).c_str());



	}
	this->m_bim.close();

	sort_data::sort((const void *)m_bimf_snpsid, m_bimf_sorted, m_num_of_snps, (DATA2SORT)D_CHARSTAR, (int)0, (int)0);

}
//=======================================================================

