## Parallel Fast AAI

A parallel version of [Fast AAI](https://github.com/cruizperez/FastAAI)

## Installation

### Pre-requisites

1. Linux (Tested with Ubuntu 22.04)
2. C++ 17 with Open MPI (Tested with gcc 11.4)
3. CMake (version 3.2)


### Building:

For a quick build, run ```build.sh``` which uses gcc and g++ to build the 
executable ```par_fastaai.x```.

For building with CMake, the following steps:

1. Clone the repository 
```
git clone https://github.com/AluruLab/ParFastAAI.git
```
2.  Create a build directory
```
mkdir ParFastAAI/build
```
3. Configure with cmake
```
cd ParFastAAI/build
cmake ..
```
4. Build with make
```
make
```
The executable ```par_fastaai.x``` will be built in the build directory


## Usage

Usage of parallel FastAAI is as follows:

    Usage: ./par_fastaai.x [OPTIONS] [path_to_query_db] [path_to_output_file]
    
    Positionals:
      path_to_query_db TEXT:FILE  Path to the Query Database
      path_to_output_file TEXT    Path to output csv file.
    
    Options:
      -h,--help                   Print this help message and exit
      -s,--separator TEXT [,]     Field Separator in the output file [Optional (default: ,)].
      -q,--query_subset TEXT:FILE Path to Query List (Should be subset of genomoes in input DB.)

Currently parallel FastAAI allows two types of usage:

1. Given a database of protein and tetramer information, compute 
   AJI (Average Jaccard Index) for all the pairs and outputs a csv file with
   AJI matrix
2. Given a database of protein and tetramer information and a query set of genomes 
   (a subset of those in the database), compute 
   AJI (Average Jaccard Index) for the query genomes against all the genomes 
   in the database and outputs a csv file with AJI matrix

data/ directory contains the example databases and the output files.
