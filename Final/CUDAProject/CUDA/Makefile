# Warnings
WFLAGS  := -Wall -Wextra -Wsign-conversion -Wsign-compare

# Optimization and architecture
OPT             := -O3
ARCH    := -march=native

# Language standard
CCSTD   := -std=c99
CXXSTD  := -std=c++11

# Linker options
LDOPT   := $(OPT)
LDFLAGS :=
BIN = "/usr/local/gcc/6.4.0/bin/gcc"
.DEFAULT_GOAL := all

.PHONY: debug
debug : OPT  := -O0 -g -fno-omit-frame-pointer -fsanitize=address
debug : LDFLAGS := -fsanitize=address
debug : ARCH :=
debug : $(EXEC)

all : sudokusolver9x9 sudokusolver16x16

sudokusolver9x9 : sudokusolver9x9.cu
	module load cuda;nvcc sudokusolver9x9.cu -o sudokusolver9x9 -lcublas -ccbin "/usr/local/gcc/6.4.0/bin/gcc"
                                
sudokusolver16x16 : sudokusolver16x16.cu
	module load cuda;nvcc sudokusolver16x16.cu -o sudokusolver16x16 -lcublas -ccbin "/usr/local/gcc/6.4.0/bin/gcc"


.PHONY: clean
clean:
	@ rm -f $(EXEC) $(OBJS)
