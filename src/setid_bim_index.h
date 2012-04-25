/*************************************************************
 *
 * Moving window project 
 * File: setid_bim_index.h	
 * Date: Dec 9, 2010
 * Author: Larissa Miropolsky
 *
 * Description:
 *   create indexing - hash table between SetID.txt file and *.Bim file
 *	 [setid_row#] = bim_row#
 **************************************************************/

#ifndef _SETID_BIM_INDEX_H        
#define _SETID_BIM_INDEX_H 
#include <fstream>  
#include <iostream> 
#include <algorithm>
#include <vector>
#include <climits>
#include <cstring>
#include <memory>

#include "DArray.h"
#include "NPsort.hpp"
#include "error_messages.h"
//===============================================================
//===============================================================


typedef struct 
{
	char snp_id[SNP_ID_SIZE];
	char letters[2];// = {NULL,NULL};
	int total_counter_per_letter[2]; //= {0,0};
	int line_counter_per_letter[2]; //= {0,0};
	int flag;
} SNP_info;


//===============================================================
//===============================================================

class Hasht   // info for column i  the columns i = 0 ... get_noc()-1
{
private:
	int *index;
	
	char* m_setidfile; 
	char* m_bimfile;

	char** m_bimf_snpsid;
	int* m_bimf_sorted;


	std::ofstream m_log;
	std::ifstream m_setid;
	std::ifstream m_bim;
	SNP_info* m_snp_sets;
	


	void upload_snpid_from_bim(int * myerror);
	void upload_snpid_from_setid_build_hash(int * myerror);

	int binsearch(const char* source); 

	void Tokenize(const std::string& str,
                      std::vector<std::string>& tokens,
                      const std::string& delimiters = " ");

	void	Get_Num_of_SNPs_in_SetID(int * myerror);


public:	
	Hasht(char * setID, char * bim, char* mwa, int * myerror);
	
	~Hasht();

	SNP_info* get_snps_sets() {return this->m_snp_sets; }

	char** m_setidf_setid;
	int* m_hash_table;
	int m_num_of_snps;//number of snps
	int m_num_of_snps_insetid;

	int m_num_of_snps_insetid_org ; // added by SLEE 

};

#endif //_SETID_BIM_INDEX_H
