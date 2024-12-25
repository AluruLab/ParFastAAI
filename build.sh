#!/bin/sh

SQLITE_DIR=./ext/sqlite/
FMT_DIR=./ext/fmt/
CLI11_DIR=./ext/CLI11/
CEREAL_DIR=./ext/cereal/
PFAAI_EXE=par_fastaai.x
GCC=gcc
GCPP=g++

SQLITE_OBJ_CMD="$GCC -pthread $SQLITE_DIR/sqlite3.c $SQLITE_DIR/shell.c -ldl -o $SQLITE_DIR/sqlite3"
echo "$SQLITE_OBJ_CMD"
$SQLITE_OBJ_CMD
#
SQLITE_OBJ_CMD="$GCC $SQLITE_DIR/sqlite3.c -c -o $SQLITE_DIR/sqlite3.o"
echo "$SQLITE_OBJ_CMD"
$SQLITE_OBJ_CMD
#
FMT_OBJ_CMD="$GCC $FMT_DIR/src/format.cc -I $FMT_DIR/include/ -c -o $FMT_DIR/src/format.o"
echo "$FMT_OBJ_CMD"
$FMT_OBJ_CMD
#
MAIN_OBJ_CMD="$GCPP  -fopenmp -I include/ -I $SQLITE_DIR -I $FMT_DIR/include/ -I $CLI11_DIR/include/ -I $CEREAL_DIR/include/ src/main.cpp -c -o src/main.o"
echo "$MAIN_OBJ_CMD"
$MAIN_OBJ_CMD
#
PAR_FAAI_BUILD="$GCPP  -fopenmp -I $SQLITE_DIR $SQLITE_DIR/sqlite3.o $FMT_DIR/src/format.o src/main.o -ldl -o $PFAAI_EXE"
echo "$PAR_FAAI_BUILD"
$PAR_FAAI_BUILD
