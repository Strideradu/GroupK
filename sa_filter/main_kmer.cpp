/*
Calculate the LCS between sequence pairs of one fasta files
Use inverse suffix array for calculating LCP array
Calculate the LCP between each pair of suffixes
*/

#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <vector>
#include <dirent.h>
#include <errno.h>
#include <ctime>
#include <sstream>
#include <set>
#include <algorithm>
#include <unordered_map>
#include "radix.h"

using namespace std;
typedef unsigned int uint;

struct fasta
{
    string title;
    string seq;
};

void print_fa(const fasta &fa)
{
    cout << fa.title << endl;
    cout << fa.seq << endl;
}

void upper_str(string &seq)
{
    transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
}

vector<fasta> read_fa(const string &file_path)
{
    string line;
    string header = ">";
    ifstream f(file_path.c_str());
    vector<fasta> fa;
    string title;
    string seq;

    if (f.is_open())
    {
        while (getline(f, line))
        {
            if (line[0] == '>')
            {
                if (!title.empty() and !seq.empty())
                {
                    fasta one_seq;
                    one_seq.title = title;
                    upper_str(seq); // convert to upper case
                    one_seq.seq = seq;
                    //cout<<title<<"\t"<<seq.length()<<endl;
                    fa.push_back(one_seq);
                }
                title = line;
                seq.clear();
            }
            else
            {
                seq += line;
            }
        }
        fasta one_seq;
        one_seq.title = title;
        upper_str(seq);
        one_seq.seq = seq;
        //cout<<title<<"\t"<<seq.length()<<endl;
        fa.push_back(one_seq);
        f.close();
    }
    else
        cout << "Unable to open file\n";
    return fa;
}

int getdir(const string &dir, vector<string> &files)
{
    DIR *dp;
    struct dirent *dirp;
    if ((dp = opendir(dir.c_str())) == NULL)
    {
        cout << "Error(" << errno << ") opening " << dir << endl;
        //cout << "Error opening " << dir << endl;
        return errno;
    }

    while ((dirp = readdir(dp)) != NULL)
    {
        if (string(dirp->d_name) != string(".") && string(dirp->d_name) != string(".."))
        {
            files.push_back(dir + string(dirp->d_name));
        }
    }
    closedir(dp);
    return 0;
}

/*
int max(int a, int b){
    return (a>b)? a: b;
}*/

int max_array(int *A, int len)
{
    int max_num = A[0];
    for (int i = 0; i < len; i++)
    {
        if (max_num < A[i])
            max_num = A[i];
    }
    return max_num;
}

int min_array(int *A, int len)
{
    int min_num = A[0];
    for (int i = 0; i < len; i++)
    {
        if (min_num > A[i])
            min_num = A[i];
    }
    return min_num;
}

vector<string> split(const string &s, char delim)
{
    stringstream ss(s);
    string item;
    vector<string> tokens;
    while (getline(ss, item, delim))
    {
        tokens.push_back(item);
    }
    return tokens;
}

/*
int LCP(unsigned char* start, unsigned char* a, unsigned char* b, int len1){
	int len = 0;
	while(a<start+len1 && *b != '\0'){
		if(*a != *b)
			break;
		len++; a++; b++;
	}
	return len;
}
*/

unsigned int *LCP_SA(unsigned int *sa, string &seq)
{
    // return LCP array
    int seq_len = seq.length();
    unsigned int *inv_sa = new unsigned int[seq_len];
    for (int i = 0; i < seq_len; i++)
        inv_sa[sa[i]] = i;

    unsigned int l = 0;
    unsigned int *LCP = new unsigned int[seq_len - 1];
    for (int i = 0; i < seq_len; i++)
    {
        unsigned int k = inv_sa[i]; //suffix T[i..]
        if (k == 0)
        {
            if (l > 0)
                l--;
            continue;
        }
        unsigned int j = sa[k - 1];
        while (seq[i + l] == seq[j + l] && seq[i + l] != '$' && seq[i + l] != 'N')
            l++;
        LCP[k - 1] = l;
        if (l > 0)
            l--;
    }
    delete[] inv_sa;
    return LCP;
}

unsigned int get_suffix_idx(unsigned int sa_val, unsigned int left, unsigned int right, vector<unsigned int> &seq_len_array)
{
    //seq_len_array: an array storing the sequence lengths
    uint mid = (left + right) / 2;
    if (mid == left)
    {
        if (sa_val >= seq_len_array[left])
            return right;
        else
            return left;
    }
    if (sa_val < seq_len_array[mid])
    {
        right = mid;
        return get_suffix_idx(sa_val, left, right, seq_len_array);
    }
    else
    {
        left = mid;
        return get_suffix_idx(sa_val, left, right, seq_len_array);
    }
}

unsigned int *get_seq_index(unsigned int *sa, unsigned int sa_len, vector<unsigned int> &seq_len_array)
{
    unsigned int *result = new unsigned int[sa_len];
    for (unsigned int i = 0; i < sa_len; i++)
    {
        unsigned int index = get_suffix_idx(sa[i], 0, seq_len_array.size() - 1, seq_len_array);
        result[i] = index;
    }
    return result;
}

char seq_alphabet[] = "ATUGCYRSWKMBDHVN";
char rev_alphabet[] = "TAACGRYSWMKVHDBN";
unordered_map<char, char> convert({{'A', 'T'}, {'T', 'A'}, {'U', 'A'}, {'C', 'G'}, {'G', 'C'}, {'Y', 'R'}, {'R', 'Y'}, {'S', 'S'}, {'W', 'W'}, {'K', 'M'}, {'M', 'K'}, {'B', 'V'}, {'D', 'H'}, {'H', 'D'}, {'V', 'B'}, {'N', 'N'}});

char complement(char n)
{
    if (convert.find(n) != convert.end())
        return convert[n];
    else
        return ' ';
}

void rev_com(string &seq)
{
    transform(seq.begin(), seq.end(), seq.begin(), complement);
    reverse(seq.begin(), seq.end());
}

void cal_LCS(uint *seq_sa, uint *LCP_array, uint *seq_index_array, uint *com_LCP_array, uint seq_len, uint seq_num, uint cutoff, uint **LCS_result, bool rev = false)
{
    /*
    seq_sa: suffix array
    LCP_array: LCP of adjacent suffixes
    seq_index_array: the original sequence index of suffix i
    com_LCP_array: the LCP of current suffix i and all other suffixes after i
    seq_len: the length of the whole sequence
    seq_num1: the number of sequence in the first fasta file
    cutoff: LCS cutoff

    LCS_result: matrix to save the results
    */
    long int lcp_count = 0;
    for (uint i = 0; i < seq_len - 2; i++)
    {
        // update com_LCS_array
        uint j = i;
        uint min_len = LCP_array[j];
        if (min_len < cutoff)
            continue;
        if (i > 0)
        {
            while (j < seq_len - 2 && min_len > LCP_array[i - 1] && min_len >= cutoff)
            {
                if (min_len > LCP_array[j])
                    min_len = LCP_array[j];
                com_LCP_array[j] = min_len; //LCP of suffix [i,j+1]
                j++;
            }
        }

        j = i;
        uint seq_idx1 = seq_index_array[i];
        while (j < seq_len - 2 && com_LCP_array[j] >= cutoff)
        {
            uint seq_idx2 = seq_index_array[j + 1];
            if (rev)
            {
                if (seq_idx1 < seq_num)
                {
                    if (seq_idx2 >= seq_num)
                    {
                        uint real_idx2 = seq_idx2 - seq_num;
                        LCS_result[seq_idx1][real_idx2]++;
                    }
                }
                else
                {
                    if (seq_idx2 < seq_num)
                    {
                        uint real_idx1 = seq_idx1 - seq_num;
                        LCS_result[real_idx1][seq_idx2]++;
                    }
                }
            }
            else
            {
                lcp_count++;
                if (seq_idx1 != seq_idx2)
                {
                    LCS_result[seq_idx1][seq_idx2]++;
                }
            }

            /*
            if(rev){
                if(seq_idx1<seq_num){
                    if(seq_idx2>=seq_num){
                        uint real_idx2 = seq_idx2 - seq_num;
                        if(seq_idx1<real_idx2)
                            LCS_result[real_idx2][seq_idx1]++;
                        else if(seq_idx1>real_idx2)
                            LCS_result[seq_idx1][real_idx2]++;
                    }
                }
                else{
                    if(seq_idx2<seq_num){
                        uint real_idx1= seq_idx1 - seq_num;
                        if(real_idx1<seq_idx2)
                            LCS_result[seq_idx2][real_idx1]++;
                        else if(real_idx1>seq_idx2)
                            LCS_result[real_idx1][seq_idx2]++;
                    }
                }
            }
            else{
                if(seq_idx1!=seq_idx2){
                    if(seq_idx1<seq_idx2)
                        LCS_result[seq_idx1][seq_idx2]++;
                    else
                        LCS_result[seq_idx2][seq_idx1]++;
                    //loc[0][seq_idx1][seq_idx2] = seq_sa[i];
                }
            }
            */

            //if(rev) loc[1][seq_idx1][seq_idx2] = 1; //The LCS is between rev_com(virus_seq) and bacteria_seq
            /*
                if(LCS_result[seq_idx1][seq_idx2]<com_LCP_array[j]){
                    LCS_result[seq_idx1][seq_idx2] = com_LCP_array[j];
                    loc[0][seq_idx1][seq_idx2] = seq_sa[i];
                    if(rev) loc[1][seq_idx1][seq_idx2] = 1; //The LCS is between rev_com(virus_seq) and bacteria_seq
                }*/
            j++;
        }
    }
    cout<<"lcp array search"<<lcp_count<<endl;
}

int main(int argc, char *argv[])
{
    clock_t start_time = clock();
    char *fa_file;
    uint cutoff = 0;
    char *out_file;
    if (argc < 7)
    {
        cout << "Usage is -i <infile> -k <kmersize> -o <outfile>\n";
        cin.get();
        exit(0);
    }
    else
    {
        for (int i = 1; i < argc; i++)
        {
            if (i + 1 < argc)
            {
                if (strcmp(argv[i], "-i") == 0)
                {
                    fa_file = argv[i + 1];
                    i++;
                }
                else if (strcmp(argv[i], "-k") == 0)
                {
                    char *cut_tmp = argv[i + 1];
                    cutoff = atoi(cut_tmp);
                    i++;
                }
                else if (strcmp(argv[i], "-o") == 0)
                {
                    out_file = argv[i + 1];
                    i++;
                }
                else
                {
                    cout << argv[i] << endl;
                    cout << "Not enough or invalid arguments, please try again.\n";
                    exit(0);
                }
            }
        }
    }

    ofstream ofile;
    ofile.open(out_file);
    string str_fa_file1(fa_file);
    vector<fasta> fa1;
    vector<uint> seq_len_array;
    uint seq_idx = 0; //The total sequence length
    string seq_whole = "";
    string rev_seq = ""; //Concatenate the reverse complement sequences of the first fasta file

    if (str_fa_file1.substr(str_fa_file1.length() - 3, 3) == string(".fa") || str_fa_file1.substr(str_fa_file1.length() - 6, 6) == string(".fasta") || str_fa_file1.substr(str_fa_file1.length() - 4) == string(".fna"))
    {
        fa1 = read_fa(str_fa_file1);
        cout << fa1.size() << endl;

        for (uint i = 0; i < fa1.size(); i++)
        {
            seq_whole += fa1[i].seq + "$";
            seq_idx += fa1[i].seq.length() + 1;
            rev_com(fa1[i].seq);
            rev_seq += fa1[i].seq + "$";
            seq_len_array.push_back(seq_idx);
        }
        rev_seq += seq_whole; //add the plus sequence at the end
    }
    else
        return 0;

    int virus_seq_len = seq_idx;
    cout << "The total length of the fa file is:\t" << seq_idx << endl;

    unsigned char *seq = new unsigned char[seq_whole.length() + 1]; //'\0' as the string end
    strcpy((char *)seq, seq_whole.c_str());
    uint seq_len = seq_whole.length();                                             // total sequence length
    unsigned int *seq_sa = Radix(seq, seq_len).build();                            // suffix array
    unsigned int *seq_index_array = get_seq_index(seq_sa, seq_len, seq_len_array); // original sequence index
    unsigned int *LCP_array = LCP_SA(seq_sa, seq_whole);                           // LCP array, length: seq_len - 1

    uint seq_num1 = fa1.size();
    unsigned int **LCS_result = new unsigned int *[seq_num1];
    for (uint i = 0; i < seq_num1; i++)
    {
        LCS_result[i] = new unsigned int[seq_num1];
        fill_n(LCS_result[i], seq_num1, 0);
    }

    // the position of the LCS region and whether it is from reverse complementary
    /*
    uint*** loc = new uint** [2];
    for(uint i=0; i<2; i++){
        loc[i] = new uint* [seq_num1];
        for(uint j=0; j<seq_num1; j++){
            loc[i][j] = new uint [seq_num1];
            fill_n(loc[i][j], seq_num1, 0);
        }
    }*/

    clock_t now_time = clock();
    double elapsed_secs = double(now_time - start_time) / CLOCKS_PER_SEC;
    cout << "Running time: " << elapsed_secs << endl;

    cout << "Beginning finding LCP:" << endl;

    // only look at LCP greater or equal to cutoff
    uint *com_LCP_array = new uint[seq_len - 1];
    fill_n(com_LCP_array, seq_len - 1, 0);
    uint min_len = LCP_array[0];
    com_LCP_array[0] = LCP_array[0];
    for (uint i = 1; i < seq_len - 1; i++)
    {
        if (min_len > LCP_array[i])
            min_len = LCP_array[i];
        com_LCP_array[i] = min_len;
        if (min_len == 0)
            break;
    }
    cal_LCS(seq_sa, LCP_array, seq_index_array, com_LCP_array, seq_len, seq_num1, cutoff, LCS_result, false);

    // calculate the memory usage
    unsigned long long int memory_usage = 0; // bytes
    cout << "The size of the sequence is:\t" << sizeof(char) * seq_len << endl;
    memory_usage += 2 * sizeof(char) * seq_len;
    cout << "The size of suffix array is:\t" << sizeof(uint) * seq_len << endl;
    cout << "The size of the LCS_result is:\t" << sizeof(uint) * seq_num1 * seq_num1 << endl;
    cout << "The size of the Location is:\t" << 2 * sizeof(uint) * seq_num1 * seq_num1 << endl;
    memory_usage += 3 * sizeof(uint) * seq_num1 * seq_num1;
    cout << "The size of suffix array, LCP array, index array, and com_LCP_array is:\t" << 4 * sizeof(uint) * seq_len << endl;
    cout << 4 * sizeof(uint) * seq_len / 1024 / 1024 << "MB\n";
    memory_usage += 4 * sizeof(uint) * seq_len;

    cout << "The total memory usage is:\t" << memory_usage << "bytes"
         << "\t" << memory_usage / 1024 / 1024 << "MB" << endl;

    delete[] seq;
    delete[] seq_sa;
    delete[] LCP_array;
    delete[] seq_index_array;
    delete[] com_LCP_array;

    now_time = clock();
    elapsed_secs = double(now_time - start_time) / CLOCKS_PER_SEC;
    cout << "Running time: " << elapsed_secs << endl;

    /*----------------------------------------------------------------------*/
    // reverse complement sequence
    seq = new unsigned char[rev_seq.length() + 1]; //'\0' as the string end
    strcpy((char *)seq, rev_seq.c_str());
    seq_len = rev_seq.length();
    seq_sa = Radix(seq, seq_len).build(); // suffix array
    cout << "Reverse suffix array finished!" << endl;
    uint idx_end = seq_len_array.back();
    uint array_size = seq_len_array.size();
    cout << seq_len << endl;
    for (uint i = 0; i < array_size; i++)
    {
        seq_len_array.push_back(seq_len_array[i] + idx_end);
    }
    cout << seq_len_array.size() << endl;
    seq_index_array = get_seq_index(seq_sa, seq_len, seq_len_array); // original sequence index
    LCP_array = LCP_SA(seq_sa, rev_seq);                             // LCP array, length: seq_len - 1

    cout << "Beginning finding LCP for reverse complementary sequences:" << endl;
    // only look at LCP greater or equal to cutoff
    com_LCP_array = new uint[seq_len - 1];
    fill_n(com_LCP_array, seq_len - 1, 0);
    min_len = LCP_array[0];
    com_LCP_array[0] = LCP_array[0];
    for (uint i = 1; i < seq_len - 1; i++)
    {
        if (min_len > LCP_array[i])
            min_len = LCP_array[i];
        com_LCP_array[i] = min_len;
        if (min_len < cutoff)
            break;
    }
    unsigned int **rev_LCS_result = new unsigned int *[seq_num1];
    for (uint i = 0; i < seq_num1; i++)
    {
        rev_LCS_result[i] = new unsigned int[seq_num1];
        fill_n(rev_LCS_result[i], seq_num1, 0);
    }
    cal_LCS(seq_sa, LCP_array, seq_index_array, com_LCP_array, seq_len, seq_num1, cutoff, rev_LCS_result, true);

    delete[] seq;
    delete[] seq_sa;
    delete[] LCP_array;
    delete[] seq_index_array;
    delete[] com_LCP_array;
    cout << "Calculating reverse LCS finished!" << endl;

    // preprocess the result
    for (uint m = 0; m < seq_num1 - 1; m++)
    {
        for (uint n = m + 1; n < seq_num1; n++)
        {
            LCS_result[m][n] += LCS_result[n][m];
            rev_LCS_result[m][n] += rev_LCS_result[n][m];
            if(rev_LCS_result[m][n]%2!=0) cout<<"unnormal value: "<<rev_LCS_result[m][n]<<endl;
            rev_LCS_result[m][n] /=2;
        }
    }
    // output the result
    uint min_LCS = 2147483647, max_LCS = 0;
    cout << seq_num1 << endl;
    for (uint m = 0; m < seq_num1 - 1; m++)
    {
        for (uint n = m + 1; n < seq_num1; n++)
        {
            if (LCS_result[m][n] > 0)
                ofile << fa1[m].title << "\t" << fa1[n].title << "\t+"
                      << "\t" << LCS_result[m][n] << endl;
            if (rev_LCS_result[m][n] > 0)
                ofile << fa1[m].title << "\t" << fa1[n].title << "\t-"
                      << "\t" << rev_LCS_result[m][n] << endl;
            /*
                string LCS_seq;
                if(loc[1][m][n]){
                    LCS_seq = rev_seq.substr(loc[0][m][n], LCS_result[m][n]); //LCS between rev_com(virus) and bacteria
                    ofile<<"-\t"<<LCS_seq<<endl;
                }
                else{
                    LCS_seq = seq_whole.substr(loc[0][m][n], LCS_result[m][n]);
                    ofile<<"+\t"<<LCS_seq<<endl;
                }*/

            if (min_LCS > LCS_result[m][n])
                min_LCS = LCS_result[m][n];
            if (max_LCS < LCS_result[m][n])
                max_LCS = LCS_result[m][n];
        }
    }
    ofile.close();
    cout << "The minimum kmer number is:\t" << min_LCS << endl;
    cout << "The maximum kmer number is:\t" << max_LCS << endl;

    for (uint i = 0; i < seq_num1; i++)
    {
        delete[] LCS_result[i];
        delete[] rev_LCS_result[i];
    }
    delete[] LCS_result;
    delete[] rev_LCS_result;

    /*
    for(uint i=0; i<2; i++){
        for(uint j=0; j<seq_num1; j++){
            delete [] loc[i][j];
        }
        delete [] loc[i];
    }
    delete [] loc;
    */

    now_time = clock();
    elapsed_secs = double(now_time - start_time) / CLOCKS_PER_SEC;
    cout << "Running time: " << elapsed_secs << endl;

    return 0;
}
