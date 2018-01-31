/*One classic algorithm problem is the LCS (longest common substring). In the past we
solve it using dynamic programming. Now we use suffix array to solve it.*/
//To implement suffix array based algorithm we need sorting. To sort C characters there
//are some tricks when using the std sort routine. We need a customized comparator.
#include "SuffixArray.h"

struct cmp
{
  bool operator() (char *a, char *b){
		return strcmp(a, b) < 0;
	}
};

int comLen(const char* a, const char* b){
	int len = 0;
	while(*a !='\0' && *b != '\0'){
		if(*a != *b)
			break;
		len++; a++; b++;
	}
	return len;
}

char* SuffixArray::GetLongestCommon(char *A, char *B)
{
	int lena = strlen(A), lenb = strlen(B);
	char **SuffixA = new char*[lena];
    char **SuffixB = new char*[lenb];

	//create suffix array for A
	int maxLen = 0;
	int maxIdx = 0;
	for(int i=0; i<lena; i++)
		SuffixA[i] = &A[i];
	for(int i=0; i<lenb; i++)
		SuffixB[i] = &B[i];

    //for(int i=0; i<lena; i++)
	//	cout << SuffixA[i] << " ";
	//cout << endl;

	sort(&SuffixA[0],&SuffixA[lena],cmp()); //There are three tricks here.  1. We should use &SuffixA[0], not SuffixA[0]. This is because
	sort(&SuffixB[0],&SuffixB[lenb],cmp()); //sort accepts pointer or iterator type, which points to the elements to be sorted. Here we
	int idxA=0, idxB=0;                     //want to sort char *, not char. If we put it as SuffixA[0], the things get sorted is the char
	while(idxA<lena && idxB<lenb){          //pointed to by SuffixA[0]...SuffixA[lena].  2. The end boundary is &SuffixA[lena], not &SuffixA[lena-1].
		if(maxLen<comLen(SuffixA[idxA],SuffixB[idxB])){  //This is so since sort use an exclusive right boundary. Similarly, if we sort an vector V,
			maxLen = comLen(SuffixA[idxA],SuffixB[idxB]); //we should use sort(V.begin(), V.end()), right?    3. We need a customized comparator as.
			maxIdx = idxA;                                //there is no default comparator for char *.
		}
		if(strcmp(SuffixA[idxA],SuffixB[idxB])==0){
			idxA++; idxB++;
		}
		else if(strcmp(SuffixA[idxA],SuffixB[idxB])>0)
			idxB++;
		else
			idxA++;
	}
	char *res = new char[maxLen+1];
	strncpy(res,SuffixA[maxIdx],maxLen);
	res[maxLen] = '\0';
	//for(int i=0; i<lena; i++) delete [] SuffixA[i];
	delete [] SuffixA;
	//for(int i=0; i<lenb; i++) delete [] SuffixB[i];
	delete [] SuffixB;
	return res;
}
