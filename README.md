## Parallel Fast AAI

A parallel version of [Fast AAI](https://github.com/cruizperez/FastAAI)

## Installation

### Pre-requisites

1. Linux (Tested with Ubuntu 22.04)
2. C++ 17 with Open MPI (Tested with gcc 11.4)
3. CMake (version 3.2, Optional)


### Build par_fastaai

For a quick build, clone the repo and run ```build.sh```. This build uses 
gcc and g++ to build the executable ```par_fastaai.x```.

For building with CMake, execute the following steps:

1. Clone the repository 
```
git clone --recurse-submodules https://github.com/AluruLab/ParFastAAI.git
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

[data/](data/) directory contains the example databases and the output files.


## Execution

Once the executable has been built,  the following for more information on all
the options that the executable accepts:
<pre><code>./par_fastaai.x --help
</code></pre>

By default the program uses all the cores in the machine. To reduce the number
of cores use the environment variable `OMP_NUM_THREADS` as described  in the
[OpenMP documentation](https://www.openmp.org/spec-html/5.0/openmpse50.html)
<pre><code> setenv OMP_NUM_THREADS 4,3,2 
</code></pre>  

## Parallel Algorithm

The parallel algorithm and the data structures used are described in the
document [Parallel Fast AAI](doc/pfaai_algorithm.pdf)

## Licensing

Our code is licensed under the Apache License 2.0 (see [`LICENSE`](LICENSE)).
