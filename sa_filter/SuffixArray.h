#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstring>

using namespace std;

class SuffixArray
{
public:
  SuffixArray(void) {}
	~SuffixArray(void) {}
	static char *GetLongestCommon(char *A, char *B);
};
