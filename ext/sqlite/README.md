# SQLite Source code

This directory contains the SQLite source code from <https://sqlite.org>.

## Public Domain

The SQLite source code is in the public domain.  See
<https://sqlite.org/copyright.html> for details. 


## The SQLite Amalgamation

The SQLite library consists of 111 files of C code (as of Version 3.37.0 - 2021-11-27) in the core with 22 additional files that implement certain commonly used extensions. Of the 133 main source files, about 75% are C code and about 25% are C header files. Most of these are "source" files in the sense that they are stored in the SQLite version control system and are edited manually in an ordinary text editor. But some of the C-language files are generated using scripts or auxiliary programs. For example, the parse.y file contains an LALR(1) grammar of the SQL language which is compiled, by the Lemon parser generator, to produce a parser contained in the file "parse.c" accompanied by token identifiers in "parse.h".

The makefiles for SQLite have an "sqlite3.c" target for building the amalgamation, to contain all C code for the core SQLite library and the FTS3, FTS5, RTREE, DBSTAT, JSON1, RBU and SESSION extensions. This file contains about 238K lines of code (or 145K if you omit blank lines and comments) and is over 8.4 megabytes in size (as of 2021-12-29).

Though the various extensions are included in the "sqlite3.c" amalgamation file, they are disabled using #ifdef statements. Activate the extensions using compile-time options like:

    -DSQLITE_ENABLE_FTS3
    -DSQLITE_ENABLE_FTS5
    -DSQLITE_ENABLE_RTREE
    -DSQLITE_ENABLE_DBSTAT_VTAB
    -DSQLITE_ENABLE_RBU
    -DSQLITE_ENABLE_SESSION 

The amalgamation contains everything you need to integrate SQLite into a larger project. Just copy the amalgamation into your source directory and compile it along with the other C code files in your project. (A more detailed discussion of the compilation process is available.) You may also want to make use of the "sqlite3.h" header file that defines the programming API for SQLite. The sqlite3.h header file is available separately. The sqlite3.h file is also contained within the amalgamation, in the first few thousand lines. So if you have a copy of sqlite3.c but cannot seem to locate sqlite3.h, you can always regenerate the sqlite3.h by copying and pasting from the amalgamation.

In addition to making SQLite easier to incorporate into other projects, the amalgamation also makes it run faster. Many compilers are able to do additional optimizations on code when it is contained with in a single translation unit such as it is in the amalgamation. We have measured performance improvements of between 5 and 10% when we use the amalgamation to compile SQLite rather than individual source files. The downside of this is that the additional optimizations often take the form of function inlining which tends to make the size of the resulting binary image larger.



