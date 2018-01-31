GroupK
====

``GroupK`` is a overlap detection tool for DNA sequnces from PacBio. 

Usage
----------
* Python 3 (You may need to modify some script if you want to use Python 2.7)
* [CMake](https://cmake.org/) 3.7 or higher 
* GCC 4.8 or higher
* [Biopython](http://biopython.org/) library for Python
* [SortedContainers](http://www.grantjenks.com/docs/sortedcontainers/index.html) library fot Python

Usage
----------
At this stage, our tool consists three part:
####Suffix array filter
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

####Modified YASS for group hits
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
####Python script


References
----------

how to cite this tool:

    Du N., Chen J., Sun Y., Improve the sensitivity of detecting long read overlaps using grouped short k-mer matches, submitted to ISMB 2018