echo gcc -pthread sqlite3.c shell.c -ldl -o sqlite3
gcc -pthread sqlite3.c shell.c -ldl -o sqlite3

echo gcc sqlite3.c -c
gcc sqlite3.c -c
