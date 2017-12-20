#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <cuda.h>
#include <curand.h>
#include <cuda_runtime.h>

#define N 9
#define n 3

void load(char *FileName, int *board) {
    FILE * a_file = fopen(FileName, "r");

    if (a_file == NULL) {
        printf("File load fail!\n"); return;
    }

    char temp;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (!fscanf(a_file, "%c\n", &temp)) {
                printf("File loading error!\n");
                return;
            }

            if (temp >= '1' && temp <= '9') {
                board[i * N + j] = (int) (temp - '0');
            } else {
                board[i * N + j] = 0;
            }
        }
    }
}

__global__ void backtracking(int *new_array, int *empty_pos, int *num_empty, int num_array, int *dev_output){
	unsigned int index = blockIdx.x * blockDim.x + threadIdx.x;
	
	while(index < num_array){
		int empty_index = 0;
		int pos, current_val,val;
		for(empty_index=0; (empty_index < num_empty[index]) && (empty_index >= 0); ){
			pos = empty_pos[index*N*N + empty_index];
			current_val = new_array[index*N*N + pos]++;
			int r_flag = 1;
			int c_flag = 1;
			int b_flag = 1;
			int row = pos/N;
			int col = pos%N;
			for(int c = 0; c < N; c++){
				val = new_array[index*N*N + row*N +c];
				if(val == current_val) r_flag = 0;
			}
			if(r_flag == 1){
				for(int r = 0; r < N; r++){
					val = new_array[index*N*N + r*N +col];
					if(val == current_val) c_flag = 0;
				}
				if(c_flag == 1){
					for(int r = n*(row/n); r < n; r++){
						for(int c = n*(col/n); c < n; c++){
							val = new_array[index*N*N + r*N + c ];
							if(val == current_val) b_flag = 0;
						}
					}
					if(b_flag == 1) empty_index = empty_index++;
				}
			}
			if(r_flag && c_flag && b_flag == 0)
			if(current_val >= 9 ){
				new_array[index*N*N + pos] = 0;
				empty_index--;
			}	
		}
		
		if(empty_index == num_empty[index]) break;
		index += blockDim.x * gridDim.x; 
	}
	for(int i= 0; i < N*N; i++){
		dev_output[i] = new_array[index*N*N + i];
	}
	
}





__global__ void Kernel1(int *pre_array, int *new_array, int num_array, int *counter, int *empty_pos, int *num_empty){
	unsigned int index = blockIdx.x * blockDim.x + threadIdx.x;
	
	while(index < num_array){
		int emptyflag = 0;
		//printf("%d",index);
		//printf("%d\n",num_array);
		for(int i = index*N*N; i < (index+1)*N*N; i++){
			if(pre_array[i] == 0){
				emptyflag = 1;
				int row = (i - index*N*N) / N;
				int col = (i - index*N*N) % N;
				
				/*To check which number could be here*/
				for(int num = 1; num <= N; num++){
					int r_flag = 1;
					int c_flag = 1;
					int b_flag = 1;
					/*check row*/
					for(int c = 0; c < N; c++){
						if(pre_array[index*N*N + row*N + c] == num){
							r_flag = 0;
						}
					}
					if(r_flag == 1){
						for(int r = 0; r < N; r++){
							if(pre_array[index*N*N + r*N + col] == num){
								c_flag = 0;
							}
						}
						if(c_flag == 1){
							for(int r = 3*(row/3); r < 3; r++){
								for(int c = 3*(col/3); c < 3; c++){
									if(pre_array[index*N*N +r*N +c] == num){
										b_flag = 0;
									}
								}
							}
							if(b_flag == 1){
								/*this number is available, copy the array*/
								int empty_index = 0;
								int next_index = atomicAdd(counter, 1);
								for(int r = 0; r < N; r++){
									for(int c = 0; c < N; c++){
										new_array[next_index*N*N + r*N + c]=pre_array[index*N*N + r*N + c];
										if(pre_array[index*N*N + r*N + c] == 0 && (r != row || c != col)){
											empty_pos[next_index*N*N + empty_index] = r*N + c;
											empty_index++;
										}
									}
								}
								new_array[next_index*N*N + row*N + col] = num;
								num_empty[next_index] = empty_index;
								
							}
						}
					}
				}
			}
			if(emptyflag == 1) break;
		}
		index += blockDim.x * gridDim.x;
	}
	
}





int main(int argc, char* argv[])
{
	FILE *fd,*ff;
	char temp;
	int  Output[81];
	char c='\n';
	if(argc < 5){
		printf("Usage: <number of threads per block> <number of blocks> <input sudoku file> <output file>\n");
		return -1;
	}
	int Blocksize = atoi(argv[1]);
	int NumBlock = atoi(argv[2]);
	
	int *pre_array,*new_array;
	int *counter;
	int *empty_pos, *num_empty;
	int *dev_output;
	
	int a = pow(2, 26);
	cudaMalloc(&pre_array, a * sizeof(int));
	cudaMalloc(&new_array, a * sizeof(int));
	cudaMalloc(&counter, sizeof(int));
	cudaMalloc(&empty_pos, a * sizeof(int));
	cudaMalloc(&num_empty, a * sizeof(int));
	cudaMalloc(&dev_output, N * N * sizeof(int));
	
	int *Input = new int[N*N];
	char* filename = argv[3];
	load(filename, Input);
	/*Read sudoku into Input array*/
/*	fd = fopen(argv[3], "r");
	if(fd == NULL){
		printf("Failed to open file: %s\n",argv[3]);
		return -1;
	}
	
	for(int i = 0; i < N; i++){
		for(int j = 0; j < N; j++){
			if(fscanf(fd,"%c\n",&temp) == -1){
				printf("Failed to read file\n");
				return -1;
			}
			if(temp < '0' || temp > '9'){
				printf("ERROR: Input file is wrong\n");
			}else{
				Input[N*i+j] = (int)(temp - '0');
			}
		}
	}
*/	
	cudaMemcpy(pre_array, Input, N*N*sizeof(int), cudaMemcpyHostToDevice);
	
	int num_array = 1;
	Kernel1<<<NumBlock, Blocksize>>>(pre_array, new_array, num_array, counter, empty_pos, num_empty);
	
	/*Loop to find all emepty position in the borad and save all new boards*/	
	for(int i = 0; i < 9; i++){
		cudaMemcpy(&num_array, counter, sizeof(int), cudaMemcpyDeviceToHost);
		printf("total boards after an iteration %d: %d\n", 2*i, num_array);
		cudaMemset(counter, 0, sizeof(int));
		Kernel1<<<NumBlock, Blocksize>>>(new_array, pre_array, num_array, counter, empty_pos, num_empty);
		
		cudaMemcpy(&num_array, counter, sizeof(int), cudaMemcpyDeviceToHost);
		printf("total boards after an iteration %d: %d\n", 2*i+1, num_array);
		cudaMemset(counter, 0, sizeof(int));
		Kernel1<<<NumBlock, Blocksize>>>(pre_array, new_array, num_array, counter, empty_pos, num_empty);
	}
	cudaMemcpy(&num_array, counter, sizeof(int), cudaMemcpyDeviceToHost);
	
	//backtracking<<<NumBlock, Blocksize>>>(new_array, empty_pos, num_empty, num_array, dev_output);
	/*print the results*/
	ff = fopen(argv[4],"w");
	if(ff == NULL){
		printf("Failed to open file: %s\n",argv[4]);
		return -1;
	}
	
	for(int i = 0; i < N; i++){
		for(int j = 0; j < N; j++){
			if(fprintf(ff,"%c\t",Input[N*i+j]) == -1){
				printf("Failed to print\n");
				return -1;
			}
		}
		if(fprintf(ff,"%c",c) == -1){
			printf("Failed to print\n");
			return -1;
		}
	}
	
	
	cudaFree(pre_array);
	cudaFree(new_array);
	cudaFree(empty_pos);
	cudaFree(num_empty);
	cudaFree(counter);
	cudaFree(dev_output);
	
	return 0;
}