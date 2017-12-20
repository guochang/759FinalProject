#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <cuda.h>
#include <curand.h>
#include <cuda_runtime.h>

#define N 9
#define n 3


__global__ void backtracking(int *new_array, int *empty_pos, int *num_empty, int num_array, int *dev_output){
	unsigned int index = blockIdx.x * blockDim.x + threadIdx.x;

	while(index < num_array){
		
		int empty_index = 0;
		int pos, current_val,val;

		for(empty_index=0; (empty_index < num_empty[index]) && (empty_index >= 0); ){
			pos = empty_pos[index*N*N + empty_index];
			new_array[index*N*N + pos]++;
			current_val = new_array[index*N*N + pos];
			//printf("%d\n",current_val);
			int r_flag = 1;
			int c_flag = 1;
			int b_flag = 1;
			int a_flag = 1;
			int row = pos/N;
			int col = pos%N;
			for(int c = 0; c < N; c++){
				if((row*N+c) != pos){
				  val = new_array[index*N*N + row*N +c];
				  if(val == current_val) r_flag = 0;
				}
			}
			if(r_flag == 1){
			    for(int r = 0; r < N; r++){
				   if((r*N+col) != pos){
				       val = new_array[index*N*N + r*N +col];
					   if(val == current_val) c_flag = 0;
			       }
			    }
				
				if(c_flag == 1){
				    int ridx = row / n;
                    int cidx = col / n;
				
				    for(int r = 0; r < n; r++){
					   for(int c = 0; c < n; c++){
						   if(( (ridx*n+r)*N + cidx*n + c) != pos){
							   val = new_array[index*N*N + (ridx*n+r)*N + cidx*n + c ];
							   if(val == current_val) b_flag = 0;
						   }
					   }
				    }
					if(b_flag == 1){
				        if(current_val > 9 ){
                           a_flag = 0;
				           }
					}
				}
			}
					
				
			
			if((r_flag == 0) || (c_flag == 0) || (b_flag == 0) || (a_flag == 0)){
			if(current_val >= 9 ){
				new_array[index*N*N + pos] = 0;
				empty_index--;
			}
			}else{
				empty_index++;
			}	
		}
		
		if(empty_index == num_empty[index]){
		
		  for(int i= 0; i < N*N; i++){
		     dev_output[i] = new_array[index*N*N + i];
	      }
		  break;
		}
		index += blockDim.x * gridDim.x; 
	}

	
}



__global__ void Kernel1(int *pre_array, int *new_array, int num_array, int *counter, int *empty_pos, int *num_empty){
	unsigned int index = blockIdx.x * blockDim.x + threadIdx.x;
	
	while(index < num_array){
		
		int emptyflag = 0;
		
		for(int i = index*N*N; i < (index * N * N) + N * N; i++){
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
							int r_b = row / n;
                            int c_b = col / n;
							for(int r = 0; r < n; r++){
								for(int c = 0; c < n; c++){
									if(pre_array[index*N*N +(r_b*n+r)*N + c_b*n + c] == num){
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
	int  *Input,*Output;
	char c='\n';
	if(argc < 5){
		printf("Usage: <number of threads per block> <number of blocks> <input sudoku file> <output file>\n");
		return -1;
	}
	int Blocksize = atoi(argv[1]);
	int NumBlock = atoi(argv[2]);
	
	int *pre_array;  /*Stores the previous version of sudoku boards */
	int *new_array;  /*Stores the new version of sudoku boards */
	int *counter;    /*Total numbers of sudoku boards*/
	int *empty_pos;  /*Stores the position of empty space*/
	int *num_empty;  /*Number of empty space*/
	int *dev_output; /*store the finished version of sudoku board*/
	
	int a = pow(2, 26);
	cudaMalloc(&pre_array, a * sizeof(int));
	cudaMalloc(&new_array, a * sizeof(int));
	cudaMalloc(&counter, sizeof(int));
	cudaMalloc(&empty_pos, a * sizeof(int));
	cudaMalloc(&num_empty, a * sizeof(int));
	cudaMalloc(&dev_output, N * N * sizeof(int));
	
	Input = (int*)malloc(N*N*sizeof(int));
    Output = (int*)malloc(N*N*sizeof(int));
    /*Read from input file*/
	fd = fopen(argv[3], "r");
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
	
    cudaMemset(counter, 0, sizeof(int));
    cudaMemset(new_array, 0, a * sizeof(int));
    cudaMemset(pre_array, 0, a * sizeof(int));
	cudaMemcpy(pre_array, Input, N*N*sizeof(int), cudaMemcpyHostToDevice);
	
    cudaEvent_t start1;
    cudaEventCreate(&start1);
    cudaEvent_t stop1;
    cudaEventCreate(&stop1);
    cudaEventRecord(start1, NULL);

	int num_array = 1;
	Kernel1<<<NumBlock, Blocksize>>>(pre_array, new_array, num_array, counter, empty_pos, num_empty);
	cudaMemcpy(&num_array, counter, sizeof(int), cudaMemcpyDeviceToHost);
	/*Loop to find all emepty position in the borad and save all new boards*/	
	for(int i = 0; i < 24; i++){
		cudaMemset(counter, 0, sizeof(int));
		if(i % 2==0){
		    Kernel1<<<NumBlock, Blocksize>>>(new_array, pre_array, num_array, counter, empty_pos, num_empty);
		}else{
		    Kernel1<<<NumBlock, Blocksize>>>(pre_array, new_array, num_array, counter, empty_pos, num_empty);
	    }
		cudaMemcpy(&num_array, counter, sizeof(int), cudaMemcpyDeviceToHost);
		printf("Number of boards created after an iteration %d: %d\n", i, num_array);
	}
	
	backtracking<<<NumBlock, Blocksize>>>(new_array, empty_pos, num_empty, num_array, dev_output);
	
	cudaEventRecord(stop1, NULL);
    cudaEventSynchronize(stop1);
    float msecTotal1 = 0.0f;
    cudaEventElapsedTime(&msecTotal1, start1, stop1);
	printf("\n***********************************");
    printf("\nThe execution time is %f ms",msecTotal1);
	printf("\n***********************************\n");
	/*print the results*/
	
	cudaMemcpy(Output, dev_output, N * N * sizeof(int), cudaMemcpyDeviceToHost);
	
	ff = fopen(argv[4],"w");
	if(ff == NULL){
		printf("Failed to open file: %s\n",argv[4]);
		return -1;
	}
	
	for(int i = 0; i < N; i++){
		for(int j = 0; j < N; j++){
			if(fprintf(ff,"%d\t",Output[N*i+j]) == -1){
				printf("Failed to print\n");
				return -1;
			}
		}
		if(fprintf(ff,"%c",c) == -1){
			printf("Failed to print\n");
			return -1;
		}
	}
	printf("The result has been written into %s\n",argv[4]);
	
	
	for (int i = 0; i < N; i++) {
        if (i % n == 0) {
            printf("-------------------------\n");
        }

        for (int m = 0; m < N; m++) {
            if (m % n == 0) {
               printf("| ");
            }
            printf("%d ", Output[i * N + m]);
        }

        printf("|\n");
    }
    printf("-------------------------\n");
	
	free(Input);
	free(Output);
	cudaFree(pre_array);
	cudaFree(new_array);
	cudaFree(empty_pos);
	cudaFree(num_empty);
	cudaFree(counter);
	cudaFree(dev_output);
	
	
	return 0;
}