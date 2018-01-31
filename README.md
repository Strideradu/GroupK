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

####Modified YASS for group hits
A prebuild excecuteble file can found under /yass/cmake-build-debug

And to use CMake to build under Linux
```
cd yass
cmake -G "Unix Makefiles"
make
```
After finished build, you can use the python script to excecute it to get output of groups of given fasta files
####Python script


References
----------

how to cite this tool:

    Du N., Chen J., Sun Y., Improve the sensitivity of detecting long read overlaps using grouped short k-mer matches, submitted to ISMB 2018