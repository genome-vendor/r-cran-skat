/*************************************************************
 *
 * NPTest Project
 * File: NPsort.hpp	
 * Date: January 1, 2011
 * Author: Larissa Miropolsky
 *
 * Description:
 *   sort char**, int*, double* array    
 **************************************************************
 * Usage: 	sort_data::sort(arr, b, size, dtype, offset, flag); //for D_CHARSTAR
 * Example:
	int b[12];
	int size = 10;
	int offset = 0;	

	char* arr[] = {"fred", "barney", "zoot", "jim", "peter","fred", "barney", "zoot", "jim", "peter"};
	DATA2SORT dtype = D_CHARSTAR;
	bool flag = 0;
	sort_data::sort(arr, b, size, dtype, offset, flag); //for D_CHARSTAR
	flag = 1;
	sort_data::sort(arr, b, size, dtype, offset, flag); //for D_CHARSTAR
	
	int arr5 [12] = {4,2,7,8,1,0,10,0,15,1,2,3};
	dtype = D_INT;
	size = 12;
	flag = 0;
	sort_data::sort(arr5, b, size, dtype, offset, flag); //for D_INT

	double arr4 [5] = {4.6,2.8,7.3,8,11};
	dtype = D_DOUBLE;
	size = 5;
	flag = 1;
	sort_data::sort(arr4, b, size, dtype, offset, flag); //for D_DOUBLE

 **************************************************************/
#ifndef _NPSORT_HPP        
#define _NPSORT_HPP 


#include <iostream>  

enum DATA2SORT
{
	D_INT,
	D_DOUBLE,
	D_CHARSTAR,
};

class sort_data  
{

public:

	static void sort(const void* a, int* b, int size, DATA2SORT type, int offset = 0, 
		int f = 0, int (*pt2Func)(const void*, const void *) = NULL);

	static void sort(const int* a, int* b, int size, int f = 0, int (*pt2Func)(const void*, const void *) = NULL,int offset = 0);
	static void sort(const double* a, int* b, int size, int f = 0, int (*pt2Func)(const void*, const void *) = NULL, int offset = 0);
	static void sort(char** a, int* b, int size, int f = 0, int (*pt2Func)(const void*, const void *) = NULL, int offset = 0);

	static void quicksort_asc(int arr[], int left, int right, const int* a, int (*pt2Func)(const void*, const void *)) ;
	static void quicksort_desc(int arr[], int left, int right, const int* a, int (*pt2Func)(const void*, const void *)) ;

	static void quicksort_asc(int arr[], int left, int right, const double* a, int (*pt2Func)(const void*, const void *)) ;
	static void quicksort_desc(int arr[], int left, int right, const double* a, int (*pt2Func)(const void*, const void *)) ;

	static void quicksort_asc(int arr[], int left, int right, char** a, int (*pt2Func)(const void*, const void *)) ;
	static void quicksort_desc(int arr[], int left, int right, char** a, int (*pt2Func)(const void*, const void *)) ;

	static int compare_strings(const void* x, const void* y);
	static int compare_ints(const void* x, const void* y);
	static int compare_doubles(const void* x, const void* y);


};

#endif //_NPSORT_HPP