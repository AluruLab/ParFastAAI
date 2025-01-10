#!/bin/sh
cmake -DCMAKE_BUILD_TYPE=Release ..
./par_fastaai_tests.x 
./par_fastaai.x data/modified_xantho_fastaai2.db test.csv 
./par_fastaai.x data/modified_xantho_fastaai2.db test.csv -q data/qsub_test_input.txt 
./par_fastaai.x data/xdb_subset1.db test.csv -r data/xdb_subset2.db 
./par_fastaai.x data/modified_xantho_fastaai2.db test.csv -q data/qsub_test_bad_input.txt 
./par_fastaai.x data/modified_xantho_fastaai2.db test.csv -r data/xdb_subset1.db 
