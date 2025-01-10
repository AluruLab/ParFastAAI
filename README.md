## Parallel Fast AAI

A parallel version of [Fast AAI](https://github.com/cruizperez/FastAAI)

## Installation

### Pre-requisites

1. Linux (Tested with Ubuntu 22.04)
1. C++ 17 with Open MPI (Tested with gcc 11.4)
1. CMake (version 3.2, Optional)

### Build par_fastaai

For a quick build, clone the repo and run `build.sh`. This build uses
gcc and g++ to build the executable `par_fastaai.x`.

For building with CMake, execute the following steps:

1. Clone the repository

```
git clone --recurse-submodules https://github.com/AluruLab/ParFastAAI.git
```

2. Create a build directory

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

The executable `par_fastaai.x` will be built in the build directory

## Usage

Usage of Parallel FastAAI is as follows:

```
Usage: ./par_fastaai.x [OPTIONS] path_to_input_db path_to_output_file

Positionals:
  path_to_input_db TEXT:FILE REQUIRED
                              Path to the Input Database
  path_to_output_file TEXT REQUIRED
                              Path to output csv file.

Options:
  -h,--help                   Print this help message and exit
  -r,--query_db TEXT:FILE     Path to the Query Database [Optional (default: Same as the Input DB)]
  -s,--separator TEXT [,]     Field Separator in the output file [Optional (default: ,)].
  -q,--query_subset TEXT:FILE Path to Query List (Should be subset of genomoes in the input DB.)
```

Currently parallel FastAAI allows three types of usage to compute 
Average Jaccard Index (AJI) for a pairs of genomes. All of them require a
SQLite database of tetramer and genome information of single copy proteins (SCP).

1. Given a database of SCPs and tetramer information, compute
   AJI for all the pairs of genomes in the database and outputs a csv file with
   AJI matrix. Output is a csv file with a square matrix, whose size is the 
   number of genomes in the database.
2. Given a database of SCPs and tetramer information and a set of query genomes
   (a strict subset of the genomes in the database), compute AJI values for
   the query genomes against all the genomes in the database.
   Output is a csv file with AJI matrix of size :
      (Number of query genomes) X (Number of total genomes)
3. Given two databases of SCP and tetramer information - a main database and 
   a query database, compute AJI (Average Jaccard Index) of the genomes in the
   query database against all the genomes in the main database.
   NOTE: NONE of the genomes in the query database should be from the input
   database. Outputs a csv file with AJI matrix of size:
      (Number of genomes in query db) X (Number of genomes in main db)

[data/](data/) directory contains the example databases, input and the output files.

## Execution

Once the executable has been built, the following for more information on all
the options that the executable accepts:

<pre><code>./par_fastaai.x --help
</code></pre>

By default the program uses all the cores in the machine. To reduce the number
of cores use the environment variable `OMP_NUM_THREADS` as described in the
[OpenMP documentation](https://www.openmp.org/spec-html/5.0/openmpse50.html)

<pre><code> setenv OMP_NUM_THREADS 4,3,2 
</code></pre>  

## Parallel Algorithm

The parallel algorithm and the data structures used are described in the
document [Parallel Fast AAI](doc/pfaai_algorithm.pdf)

## Licensing

Our code is licensed under the Apache License 2.0 (see [`LICENSE`](LICENSE)).
