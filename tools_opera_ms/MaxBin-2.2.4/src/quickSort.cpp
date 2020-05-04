#include "iostream"
#include "quickSort.h"

quickSort::quickSort()
{
	arr = NULL;
	rank = NULL;
	num = 0;
	inputchar = NULL;
}

quickSort::quickSort(double *input, int input_num)
{
	arr = input;
	rank = NULL;
	inputchar = NULL;
	num = input_num;
}

quickSort::quickSort(double *input, int *input_rank, int input_num)
{
	arr = input;
	rank = input_rank;
	inputchar = NULL;
	num = input_num;
}

quickSort::quickSort(double *input, char **input_char, int input_num)
{
	arr = input;
	rank = NULL;
	inputchar = input_char;
	num = input_num;
}

void quickSort::input(double *input, int input_num)
{
	arr = input;
	rank = NULL;
	inputchar = NULL;
	num = input_num;
}

void quickSort::input(double *input, int *input_rank, int input_num)
{
	arr = input;
	rank = input_rank;
	inputchar = NULL;
	num = input_num;
}

void quickSort::input(double *input, char **input_char, int input_num)
{
	arr = input;
	rank = NULL;
	inputchar = input_char;
	num = input_num;
}

void quickSort::sort_all()
{
	sort(0, num - 1);
}

void quickSort::sort(int low, int high)
{
	int i, j;
	if (low < high)
	{
		i = low;
		j = high + 1;
		pivot = arr[i];
		while (i < j)
		{
			i++;
			while (arr[i] > pivot && i < num)
			{
				i++;
			}
			j--;
			while (arr[j] < pivot && j > 0)
			{
				j--;
			}
			if (i < j)
			{
				swap(&arr[i], &arr[j]);
				if (rank != NULL)
				{
					swap(&rank[i], &rank[j]);
				}
				if (inputchar != NULL)
				{
					swap(&(inputchar[i]), &(inputchar[j]));
				}
			}
		}
		swap(&arr[low], &arr[j]);
		if (rank != NULL)
		{
			swap(&rank[low], &rank[j]);
		}
		if (inputchar != NULL)
		{
			swap(&(inputchar[low]), &(inputchar[j]));
		}

		sort(low, j - 1);
		sort(j + 1, high);
	}
}

void quickSort::swap(double *m, double *n)
{
	fTemp = *m;
	*m = *n;
	*n = fTemp;
}

void quickSort::swap(int *m, int *n)
{
	iTemp = *m;
	*m = *n;
	*n = iTemp;
}

void quickSort::swap(char **m, char **n)
{
	cTemp = *m;
	*m = *n;
	*n = cTemp;
}

