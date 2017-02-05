#include<omp.h>
#include<stdio.h>
#include<inttypes.h>
#include<stdint.h>
#include<string.h>
#include<stdlib.h>
#define MAX_THREADS 12
int nthr;
typedef struct node			//This is a structure that represents a node of the Cartesian tree and suffix tree
{
	unsigned long firstChar;
	unsigned long parent;
	unsigned long value;		//value=how many characters after the first character
} node;
typedef struct				//Structure representing the suffix array
{
	unsigned char* A;		//Main string
	unsigned* SArr;			//Suffix array
	unsigned length;		//Length of string
}suffix_struct;
#define BYTE_SWAP(x) ((x >> 56) | \
  (((x >> 48) & 0xff) << 8) | \
  (((x >> 40) & 0xff) << 16) | \
  (((x >> 32) & 0xff) << 24) | \
  (((x >> 24) & 0xff) << 32) | \
  (((x >> 16) & 0xff) << 40) | \
  (((x >> 8) & 0xff) << 48) | \
  (((x) & 0xff) << 56))			//Function to swap the bytes
#define FTAB_SIZE 256			//Size of frequency array, because there are 256 characters
#ifndef REC_SORT_LIMIT			//Limit of size for selection sort, to keep the time from exceeding by the order of n^2
#define REC_SORT_LIMIT 12
#endif
__thread unsigned long long value_array[REC_SORT_LIMIT+1];	//Array for values of the suffixes
static inline int compSuffixes(unsigned n, unsigned char* arr, unsigned s1, unsigned s2, unsigned depth)
{
	//Compare and return 0 if equal, 1 if greater and -1 if less
	if (s2 + depth >= n) return 1;
	else if (s1 + depth >= n) return -1;
	else if (arr[s1+depth] == arr[s2+depth]) return 0;
	else return (arr[s1+depth] > arr[s2+depth]) ? 1 : -1;
}
static inline int compOffsetSuffixes(unsigned n, unsigned char* arr, unsigned s0, unsigned s1, unsigned depth)
{
	//Compare and return 0 if equal, 1 if greater and -1 if less
	int index0 = s0 + depth >= n ? (s0 + depth - n) : s0 + depth;  
	int index1 = s1 + depth >= n ? (s1 + depth - n) : s1 + depth;
	unsigned long long v0 = *((unsigned long long*) (arr + index0));
	unsigned long long v1 = *((unsigned long long*) (arr + index1));
	v0 = BYTE_SWAP(v0);
	v1 = BYTE_SWAP(v1);
	if (v0 > v1) return 1;
	if (v0 == v1) return 0;
	return -1;
}
static inline unsigned long long getOffsetValue(unsigned n , unsigned char* arr, unsigned pre_index, unsigned depth)
{
	//Give the value at an offset of depth
	int index = pre_index + depth >= n ? (pre_index + depth - n) : pre_index + depth;  
	unsigned long long value = *((unsigned long long*) (arr + index));
	value = BYTE_SWAP(value);
	return value;
}
static inline int compWithFirst(unsigned n, unsigned char* arr, unsigned s0, unsigned long long v1, unsigned depth)
{
	//Compare and return 0 if equal, 1 if greater and -1 if less
	int index0 = s0 + depth >= n ? (s0 + depth - n) : s0 + depth;  
	unsigned long long v0 = *((unsigned long long*) (arr + index0));
	v0 = BYTE_SWAP(v0);
	if (v0 > v1) return 1;
	if (v0 == v1) return 0;
	return -1;
}
static inline int compWithSecond(unsigned n, unsigned char* arr, unsigned long long v0, unsigned s1, unsigned depth)
{
	//Compare and return 0 if equal, 1 if greater and -1 if less
	int index1 = s1 + depth >= n ? (s1 + depth - n) : s1 + depth;
	unsigned long long v1 = *((unsigned long long*) (arr + index1));
	v1 = BYTE_SWAP(v1);
	if (v0 > v1) return 1;
	if (v0 == v1) return 0;
	return -1;
}
void sub_sort(suffix_struct* ss, int start, int end, unsigned depth)
{
	unsigned* suffix_array = ss->SArr;
	unsigned char* arr = ss->A;
	unsigned n = ss->length;
	int  init_start = start, init_end = end;
	if (start >= end)		//Because the size of the array is 0
	{
		return;
	}
	else
	{
		if (end - start < REC_SORT_LIMIT)
		{
			if (select_sort(ss, start, end, depth))	//Send the portion to get selection sort if it is small
				return;
		}
		unsigned middle_index = start + (end - start) / 2;	//3 way quicksort begins if the size is big
		unsigned equal_start = start;
		unsigned equal_end   = end;
		unsigned long long median_value = 127;
		int corrected_middle_index = suffix_array[middle_index] + depth >= n ? (suffix_array[middle_index] + depth - n) : suffix_array[middle_index] + depth;			//if this function has a depth greater than 1, it is accounted for in the middle index
		unsigned long long middle_value = *((unsigned long long*) (arr + corrected_middle_index));
		middle_value = BYTE_SWAP(middle_value);
		median_value = middle_value;
		#define swap(x, y) {\
		unsigned tmp = suffix_array[x];\
		suffix_array[x] = suffix_array[y];\
		suffix_array[y] = tmp;\
		}
		while (start < end)
		{
			unsigned long long start_value = getOffsetValue(n, arr, suffix_array[start], depth); 
			if (start_value < median_value)
			{
				start++;		//It goes into the "less than" part
			}
			else if (start_value == median_value)	//Goes into "equal to" part
			{
				swap(equal_start, start);
				start++;
				equal_start++;
			}
			else			//Greater
			{
				int comp = compWithSecond(n, arr, median_value, suffix_array[end], depth);
				//compare the end of the array with the median 
				while (comp != 1)
				{
					if (comp == 0)		//End is equal
					{
						swap(equal_end, end);
						end--;
						equal_end--;
					}
					else		//End is greater
					{
						end--;
					}
					if (start >= end)
						break;
					comp = compWithSecond(n, arr, median_value, suffix_array[end], depth);
					//Again compare the median value with the end value in the array
				}
				if (start < end)		//Swap until start goes beyond end
				{
					swap(start,end);
					start++;
					end--;
				}
			}
		}
		int comp = compWithFirst(n, arr, suffix_array[start], median_value, depth);
		//Compare starting array value with median
		if (comp == 0)		//Last equal value
		{
			swap(equal_start, start);
			equal_start++;
			start++;
		}
		else if (comp == 1)	//Last greater value
		{
			end--;
		}
		else			//Last less value
		{
			start++;
		}
		//Moving elt equal to median into position
		unsigned i,j;
		if (equal_start != start)		//Start swapping
			for (i = init_start, j = start - 1; i < equal_start; ++i, --j)
				swap(i, j);
		unsigned k, l;
		if (end != equal_end)			//Start swapping
			for (k = start, l = init_end; l > equal_end; ++k, --l)
				swap(k, l);
		unsigned left_equal_count = equal_start - init_start;
		unsigned right_equal_count = init_end - equal_end;
		sub_sort(ss, init_start, start - left_equal_count - 1, depth);		//Sort "less than" elements 
		sub_sort(ss, start - left_equal_count, start + right_equal_count - 1, depth + 8);	//Sort "equal to" elements
		sub_sort(ss, start + right_equal_count, init_end, depth);		//Sort "greater than" elements
	}
}
void createSuffixArray(suffix_struct* ss)
{
	unsigned n = ss->length;
	unsigned char* arr = ss->A;
	unsigned *suffix_array = ss->SArr;
	unsigned ftab[FTAB_SIZE+1] = {0};	//Array of frequency of different characters
	unsigned i;
	for (i = 0; i < n; ++i)			//Calculating frequencies
	{
		ftab[arr[i]] += 1; 
	}
	int threadWork[MAX_THREADS+1];	//Array to store starting points of indices for different threads
	threadWork[0] = 1;
	unsigned distWork =(n%nthr==0?n/nthr:n/nthr+1);	//Calculate optimal work to be done by a thread
	int work = ftab[0];
	int curr_thr = 0;
	for (i = 1; i < FTAB_SIZE + 1; i++)	//Loop calculating the starting points
	{
		work += ftab[i];
		ftab[i] += ftab[i-1];		//Cumulative frequencies to find the occurrences easily
		if (work >= distWork)//If work assigned to the thread is at least equal to the optimal work to be done
		{
			threadWork[curr_thr+1]=i;
			work = 0;
			curr_thr++;
		}
	}
	threadWork[curr_thr+1] = FTAB_SIZE + 1;	//The last index that is possible
	for (i = 0; i < n; i++)			//Bucket sort according to frequency
	{
		unsigned byte = arr[i];
		suffix_array[ftab[byte] - 1] = i;
		ftab[byte]--;
	}
	#pragma omp parallel num_threads(nthr)		//Divide buckets among threads
	{
		int i=omp_get_thread_num(),start=threadWork[i],end=threadWork[i+1],d=1;
		for (i=start;i<end;i++)
		{
			if (ftab[i]-ftab[i-1]>0)	//If there is an occurrence of that character
				sub_sort(ss,ftab[i-1],ftab[i]-1,d);
		}
	}
}
int select_sort(suffix_struct* ss, int start, int end, unsigned depth)
{
	unsigned i, j;
	#pragma omp parallel for num_threads(nthr)
	for (i = start; i <= end; ++i)			//Fill in the values of characters
	{
		value_array[i - start] = BYTE_SWAP((*((unsigned long long*) (ss->A + ss->SArr[i]))));
	}
	for (i = start; i < end; ++i)			//Simple selection sort
	{
		int min_index = i;
		unsigned long long min_value = value_array[i - start];
		for (j = i+1; j <= end; j++)
		{
			unsigned long long current_value = value_array[j - start];
			if (current_value < min_value)
			{
				min_index = j;
				min_value = current_value;
			}
			else if (current_value == min_value)
			{
				return 0;
			}
		}
		int tmp = ss->SArr[i];
		value_array[min_index - start] = value_array[i - start];
		ss->SArr[i] = ss->SArr[min_index];
		ss->SArr[min_index] = tmp;
	}
	return 1;			//If sort is successful
}
inline unsigned long getRoot(node* nodes, ulong i)		//Get the root of the tree
{
	unsigned long root = nodes[i].parent;
	while(root!= 0 && nodes[nodes[root].parent].value == nodes[root].value)
		root = nodes[root].parent;
	return root;
}
void merge(node* N, unsigned long left, unsigned long right)		//Merge two cartesian subtrees with roots left and right
{
	unsigned long head;
	if (N[left].value > N[right].value)
	{
		head = left;
		left = N[left].parent;
	}
	else
	{
		head = right;
		right= N[right].parent;
	}
	while(1)
	{
		if (left == 0)
		{
			N[head].parent = right;
			break;
		}
		if (right == 0)
		{
			N[head].parent = left;
			break;
		}
		if (N[left].value > N[right].value)
		{
			N[head].parent = left;
			left = N[left].parent;
		}
		else
		{
			N[head].parent = right;
			right = N[right].parent;
		}
		head = N[head].parent;
	}
}
void cartesianTree(node* Nodes, unsigned long n)	//Function to construct a Cartesian Tree
{
	unsigned long tree_size,middle,start,step;
	for(step = 2; step < n * 2; step *= 2)		//Subtrees created for particular step sizes and then merged
	{
		#pragma omp parallel for num_threads(nthr) private(tree_size,middle)		//Parallelized for a particular step size
		for(start = 0; start < n; start += step)
		{
			if (n - start >= step)
			{
				tree_size = step;
			}
			else
			{
				tree_size = n - start;
			}
			if (tree_size == 2)
			{
				if (Nodes[start].value > Nodes[start + 1].value)	
				{
					Nodes[start].parent = start + 1;
				}
				else
				{
					Nodes[start + 1].parent = start;
				}
			}
			else if (tree_size > step / 2)	
			{
				middle = start + (step / 2) - 1;
				merge(Nodes, middle, middle + 1);
			}
		}
	}
}
node* suffixArrayToTree(unsigned *SA, int* LCP,int n,char *str)
{
	long m = 2*n,i;		//m is number of nodes
	node* nodes = (node *)malloc(m*sizeof(node));
	#pragma omp parallel for num_threads(nthr)
	for(i=1; i<n; i++)		//Assign values to nodes using LCP Array
	{ 
		nodes[2*i].value = LCP[i-1];
		nodes[2*i+1].value = n-SA[i]+1; //length of string including 1 past end
		nodes[2*i].parent = i;
		nodes[2*i+1].parent = i;
	}
	nodes[0].value = 0;
	nodes[1].value = n-SA[0]+1;
	nodes[0].parent = 0;
	nodes[1].parent = 0;
	free(LCP);
	cartesianTree(nodes,1);
	#pragma omp parallel for num_threads(nthr)
	for(i = 1; i < m; i++)
		nodes[i].parent = getRoot(nodes, i);
	for(i=1;i<m;i++)
	{
		unsigned long start=(i&1 == 1)?(n-nodes[i].value+1):(n-nodes[i-1].value+1);
		unsigned long offset=start+nodes[nodes[i].parent].value;
		nodes[i].firstChar = str[offset];
	}
	return nodes;
}
void arrayToLCP(char *str,unsigned *sa,int *lcp,int n)		//Function to implement Kasai's algorithm
{
	int k=0,i;
	int *rank=(int *)calloc(n,sizeof(int));			//Array for ranks of suffixes
	#pragma omp parallel for num_threads(nthr)
	for(i=0;i<n; i++)
		rank[sa[i]]=i;
	for(i=0;i<n;i++,k?k--:0)
	{
		if(rank[i]==n-1)
		{
			k=0;
			continue;
		}
		int j=sa[rank[i]+1];
		while(i+k<n && j+k<n && str[i+k]==str[j+k])	//Compare characters until an unequal character is found
			k++;
		lcp[rank[i]]=k;
	}
}
int main(int argc, char *argv[])
{
	int size=atoi(argv[1]),i;
	nthr=atoi(argv[2]);
	suffix_struct s;
	srandom(time(NULL));
	int *lcp=(int *)calloc(size,sizeof(int));
	s.A=(char *)malloc((size+1)*sizeof(char));
	#pragma omp parallel for num_threads(nthr)
	for(i=0;i<size;i++)				//Assign random characters
		s.A[i]=(unsigned char)(rand()%26+97);
	s.A[size]='\0';
	s.SArr=(unsigned *)malloc(size*sizeof(unsigned));
	s.length=size;
	double start=omp_get_wtime();
	createSuffixArray(&s);
	/*for(i=0;i<size;i++)
		printf("%d\n",s.SArr[i]);*/
	arrayToLCP(s.A,s.SArr,lcp,size);
	/*for(i=0;i<size;i++)
		printf("%d\n",lcp[i]);*/
	node *nodes=suffixArrayToTree(s.SArr,lcp,size,s.A);
	double end=omp_get_wtime();
	printf("%f\n",(end-start));
	free(s.A);
	free(s.SArr);
	free(nodes);
	return 0;
}
