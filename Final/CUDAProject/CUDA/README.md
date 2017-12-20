9x9 and 16x16 GPU based sudokusolver
----------------------------------------------------

##Compilation 

The whole solver consist of sudokusolver9x9.cu, sudokusolver16x16.cu, Makefile,
README.md and Report about the code and the GPU acceleration analysis. 

To compile 9x9 sudokusolver
'make sudokusolver9x9'
or compile with nvcc 'nvcc sudokusolver9x9.cu -o sudokusolver9x9 -lcublas -ccbin "/usr/local/gcc/6.4.0/bin/gcc"'

To compile 16x16 sudokusolver
'make sudokusolver16x16'
or compile with nvcc 'nvcc sudokusolver16x16.cu -o sudokusolver16x16 -lcublas -ccbin "/usr/local/gcc/6.4.0/bin/gcc"'

##Running

In submit.sh, change the command to: <sudokusolvername> <number of threads per block> <number of blocks> <input sudoku file directory>/<input sudoku file> <output file>

For example,
./sudokusolver9x9 16 16 9x9puzzle/puzzle9_01.txt puzzle9.txt

or

./sudokusolver16x16 16 16 16x16puzzle/puzzle16_01.txt puzzle16.txt

Then run sbatch submit.sh, the output file is in current directory.

##Notices

1. The result of the sudoku will not only show in the terminal, but also store the result into an output file. 

2. There is a non-bug bug in this program. If we input a sudoku puzzles whose empty spaces are less than 24 (this number is depend on the source code and it is changeable), this program cannot solve this puzzle.

3. This program uses cublas in the verify sudoku part. And this part will be called automatically after solving soduku part. 
