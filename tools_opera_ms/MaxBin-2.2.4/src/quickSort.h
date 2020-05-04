#ifndef __QUICKSORT_H__
#define __QUICKSORT_H__

class quickSort
{
	public:
		quickSort();
		quickSort(double *input, int input_num);
		quickSort(double *input, int *input_rank, int input_num);
		quickSort(double *input, char **input_char, int input_num);
		void input(double *input, int input_num);
		void input(double *input, int *input_rank, int input_num);
		void input(double *input, char **input_char, int input_num);
		void sort(int low, int high);
		void sort_all();

	private:
		// Variables
		double *arr;
		int *rank;
		char **inputchar;
		int num;
		int iTemp;
		double fTemp;
		char *cTemp;
		double pivot;
		// Function
		int qPartition(int low, int high);
		void swap(double *m, double *n);
		void swap(char **m, char **n);
		void swap(int *m, int *n);
};


#endif

