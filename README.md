GroupK
====

``GroupK`` is an overlap detection tool for DNA sequnces from PacBio. 
Here is the old code for the experiment. The C++ implementation will be released soon

Dependecies
----------

* Python 3 (You may need to modify some script if you want to use Python 2.7)
* [CMake](https://cmake.org/) 3.7 or higher 
* GCC 4.8 or higher
* [Biopython](http://biopython.org/) library for Python
* [SortedContainers](http://www.grantjenks.com/docs/sortedcontainers/index.html) library fot Python

Usage
----------
At this stage, our tool consists three part:

#### Suffix array filter

A prebuild excecuteble file using GCC 4.8.2 can found under /sa_filter/

Use GCC 4.8 or higher to build
```
cd sa_filter
g++ -std=c++11 *.cpp -o filter 
```
To excecute (although this will be included in the script)
```
./filter -i <infile> -k <kmersize> -o <outfile>
```
The program will do all against all to compare shared kmer

#### Modified YASS for group hits

A prebuild excecuteble file using GCC 4.8.2 can found under /yass/

And to use CMake to build under Linux
```
cd yass
cmake -G "Unix Makefiles"
make
```
After finished build, you can use the python script to excecute it to get output of groups of given fasta files

The output format is (seperated by tab)
```
query_id, query length, target_id, target_length, groups(x, diagonal, length)
```

#### Python script

To run GroupK, after build the two program above
```
cd script
python GroupK.py [-h] [--k1 K1] [--threshold THRESHOLD] [--k2 K2]
                 [--accuracy ACCURACY] [--gap GAP] [--chain CHAIN]
                 [--groupc GROUPC] [--idbase IDBASE] [--ratio RATIO]
                 [--size SIZE] [--large LARGE]
                 input output
                 
positional arguments:
    input                 path of input fasta file
    output                path of output file

optional arguments:
    -h, --help            show this help message and exit
    --k1 K1               kmer size for filtration (default: 15)
    --threshold THRESHOLD
                        count threshold for shared k1 (default: 2)
    --k2 K2               kmer size for group (default: 9)
    --accuracy ACCURACY   accuracy of the reads (default: 0.85)
    --gap GAP             gap rate (default: 0.12)
    --chain CHAIN         number of kmer in the chain to report (default: 3)
    --groupc GROUPC       the coefficeint c in paper to control chaining
                        threshold, must be non-zero float (default: 4.0)
    --idbase IDBASE       if the number of identical base meet this threshold
                        the final chain threshold will be release by the --large times (default: 400)
    --ratio RATIO         ratio of two overlap region threshold (default: 0.5)
    --size SIZE           group size threshold (default: 12)
    --large LARGE         release threhold for large overlap (default: 2)

```
Right now, the main bottle neck is the yass part.


References
----------

how to cite this tool:

    Du N., Chen J., Sun Y., Improving the sensitivity of detecting long read overlaps using grouped short k-mer matches, submitted to ECCB 2018
