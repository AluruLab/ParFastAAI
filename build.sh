#!/bin/sh

SQLITE_DIR=./ext/sqlite/

SQLITE_OBJ_CMD="gcc -pthread $SQLITE_DIR/sqlite3.c $SQLITE_DIR/shell.c -ldl -o $SQLITE_DIR/sqlite3"
echo "$SQLITE_OBJ_CMD"
$SQLITE_OBJ_CMD
#
SQLITE_OBJ_CMD="gcc $SQLITE_DIR/sqlite3.c -c -o $SQLITE_DIR/sqlite3.o"
echo "$SQLITE_OBJ_CMD"
$SQLITE_OBJ_CMD
#
MAIN_OBJ_CMD="g++-11  -fopenmp -I $SQLITE_DIR src/main.cpp -c -o src/main.o"
echo "$MAIN_OBJ_CMD"
$MAIN_OBJ_CMD
#
PAR_FAAI_BUILD="g++-11  -fopenmp -I $SQLITE_DIR $SQLITE_DIR/sqlite3.o src/main.o -ldl -o main.exe"
echo "$PAR_FAAI_BUILD"
$PAR_FAAI_BUILD
