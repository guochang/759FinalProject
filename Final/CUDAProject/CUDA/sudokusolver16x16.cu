#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <cuda.h>
#include <curand.h>
#include <cuda_runtime.h>
#include "cublas_v2.h"

#define N 16
#define n 4

void check(int *puzzle_sol){
	cudaError_t cudaStat;
	cublasStatus_t stat;
	cublasHandle_t handle;

	int i, j;
	float* A;//store the puzzle solution
	float* A_sub;//
	float* x;//multiply with the matrix to get the sum
	float* y;//store the sum->difference
	float checksum[16] = { 136, 136, 136, 136, 136, 136, 136, 136, 136, 136, 136, 136, 136, 136, 136, 136 };

	float zero;

	A = (float*)malloc(N*N*sizeof(float));
	A_sub = (float*)malloc(N*N*sizeof(float));
	x = (float*)malloc(N*sizeof(float));
	y = (float*)malloc(N*sizeof(float));

	for (i = 0; i < N; i++){
		for (j = 0; j < N; j++){
			A[i*N + j] = (float)puzzle_sol[i*N + j];

			int sub_i, sub_j;
			sub_i = (i / 4) * 4 + (j / 4);
			sub_j = (i % 4 * 4) + (j % 4);
			A_sub[sub_i*N + sub_j] = (float)puzzle_sol[sub_i*N + sub_j];
		}
	}

	for (i = 0; i < N; i++)
		x[i] = 1.0f;

	for (i = 0; i < N; i++)
		y[i] = 0.0f;


	float* d_A;
	float* d_A_sub;
	float* d_x;
	float* d_y;
	float* d_checksum;

	cudaStat = cudaMalloc((void**)&d_A, N*N*sizeof(float));
	cudaStat = cudaMalloc((void**)&d_A_sub, N*N*sizeof(float));
	cudaStat = cudaMalloc((void**)&d_x, N*sizeof(float));
	cudaStat = cudaMalloc((void**)&d_y, N*sizeof(float));
	cudaStat = cudaMalloc((void**)&d_checksum, N*sizeof(float));


	stat = cublasCreate(&handle);
	stat = cublasSetMatrix(N, N, sizeof(float), A, N, d_A, N);
	stat = cublasSetMatrix(N, N, sizeof(float), A_sub, N, d_A_sub, N);
	stat = cublasSetVector(N, sizeof(float), x, 1, d_x, 1);
	stat = cublasSetVector(N, sizeof(float), y, 1, d_y, 1);
	stat = cublasSetVector(N, sizeof(float), checksum, 1, d_checksum, 1);

	float alpha = 1.0f;
	float beta = 0.0f;
	float minus = -1.0f;
	//d_y = al*d_a *d_x + bet *d_y
	stat = cublasSgemv(handle, CUBLAS_OP_N, N, N, &alpha, d_A, N, d_x, 1, &beta, d_y, 1);//d_y will be sum of the row, 45

	stat = cublasGetVector(N, sizeof(float), d_y, 1, y, 1);
	/*printf("Sum of each row...\n");
	for (i = 0; i < N; i++){
		printf("%f\n", y[i]);
	}*/
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//d_y = d_y + minus*d_checksum
	stat = cublasSaxpy(handle, N, &minus, d_checksum, 1, d_y, 1);//d_y-checksum, 0
	stat = cublasGetVector(N, sizeof(float), d_y, 1, y, 1);
	/*printf("Difference between correct sum...\n");
	for (i = 0; i < N; i++){
		printf("%f\n", y[i]);
	}*/
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//zero=|d_y[0]|+...+|d_y[9]|
	stat = cublasSasum(handle, N, d_y, 1, &zero);

	if (zero == 0){
		printf("Row is correct!\n");
	}
	else{
		printf("Row is incorrect\n");
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//d_y = al*d_a *d_x + bet *d_y
	stat = cublasSgemv(handle, CUBLAS_OP_T, N, N, &alpha, d_A, N, d_x, 1, &beta, d_y, 1);//d_y will be sum of the column, 45

	stat = cublasGetVector(N, sizeof(float), d_y, 1, y, 1);
	/*printf("Sum of each column...\n");
	for (i = 0; i < N; i++){
		printf("%f\n", y[i]);
	}*/
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//d_y = d_y + minus*d_checksum
	stat = cublasSaxpy(handle, N, &minus, d_checksum, 1, d_y, 1);//d_y-checksum, 0

	stat = cublasGetVector(N, sizeof(float), d_y, 1, y, 1);
	/*printf("Difference between correct sum...\n");
	for (i = 0; i < N; i++){
		printf("%f\n", y[i]);
	}*/
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//zero=|d_y[0]|+...+|d_y[9]|
	stat = cublasSasum(handle, N, d_y, 1, &zero);

	if (zero == 0){
		printf("Column is correct!\n");
	}
	else{
		printf("Column is incorrect\n");
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//d_y = al*d_a *d_x + bet *d_y
	stat = cublasSgemv(handle, CUBLAS_OP_T, N, N, &alpha, d_A_sub, N, d_x, 1, &beta, d_y, 1);//d_y will be sum of the sub box, 45

	stat = cublasGetVector(N, sizeof(float), d_y, 1, y, 1);
	/*printf("Sum of each sub box...\n");
	for (i = 0; i < N; i++){
		printf("%f\n", y[i]);
	}*/
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//d_y = d_y + minus*d_checksum
	stat = cublasSaxpy(handle, N, &minus, d_checksum, 1, d_y, 1);//d_y-checksum, 0
	stat = cublasGetVector(N, sizeof(float), d_y, 1, y, 1);
	/*printf("Difference between correct sum...\n");
	for (i = 0; i < N; i++){
		printf("%f\n", y[i]);
	}*/
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//zero=|d_y[0]|+...+|d_y[9]|
	stat = cublasSasum(handle, N, d_y, 1, &zero);

	if (zero == 0){
		printf("Sub box is correct!\n");
		printf("\nThe output 16x16 sudoku result is correct.\n");
	}
	else{
		printf("Sub box is incorrect\n");
	}

	cudaFree(d_A);
	cudaFree(d_A_sub);
	cudaFree(d_x);
	cudaFree(d_y);
	cudaFree(d_checksum);
	cublasDestroy(handle);
	free(A);
	free(x);
	free(y);
	return ;
}

/* 
 * Description:
 *    This function uses backtracking to fill all empty spaces which is not complete with BFS.
 *
 */
__global__ void backtracking(int *new_array, int *empty_pos, int *num_empty, int num_array, int *dev_output){
	unsigned int index = blockIdx.x * blockDim.x + threadIdx.x;

	while(index < num_array){
		
		int empty_index = 0;
		int pos, current_val,val;

		for(empty_index=0; (empty_index < num_empty[index]) && (empty_index >= 0); ){
			/*Get the empty space's position*/
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
			/*check row*/
			for(int c = 0; c < N; c++){
				if((row*N+c) != pos){
				  val = new_array[index*N*N + row*N +c];
				  if(val == current_val) r_flag = 0;
				}
			}
			if(r_flag == 1){
				/*check column*/
			    for(int r = 0; r < N; r++){
				   if((r*N+col) != pos){
				       val = new_array[index*N*N + r*N +col];
					   if(val == current_val) c_flag = 0;
			       }
			    }
				
				if(c_flag == 1){
					/*check box*/
				    int r_b = row / n;
                    int c_b = col / n;
				
				    for(int r = 0; r < n; r++){
					   for(int c = 0; c < n; c++){
						   if(( (r_b*n+r)*N + c_b*n + c) != pos){
							   val = new_array[index*N*N + (r_b*n+r)*N + c_b*n + c ];
							   if(val == current_val) b_flag = 0;
						   }
					   }
				    }
					if(b_flag == 1){
						/*check the current value*/
				        if(current_val > 16 ){
                           a_flag = 0;
				           }
					}
				}
			}
				
			if((r_flag == 0) || (c_flag == 0) || (b_flag == 0) || (a_flag == 0)){
			if(current_val >= 16 ){
				/*backtrack to previous attempt*/
				new_array[index*N*N + pos] = 0;
				empty_index--;
			}
			}else{
				empty_index++;
			}	
		}
		
		if(empty_index == num_empty[index]){
		/*copy the result to output array*/
		  for(int i= 0; i < N*N; i++){
		     dev_output[i] = new_array[index*N*N + i];
	      }
		  break;
		}
		index += blockDim.x * gridDim.x; 
	}

	
}

/* 
 * Description:
 *    This function is used to find all possible next ruslts by using breadth
 *  first search.
 *
 */
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
						/*check column*/
						for(int r = 0; r < N; r++){
							if(pre_array[index*N*N + r*N + col] == num){
								c_flag = 0;
							}
						}
						if(c_flag == 1){
							/*check box*/
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
											/*find the position of empty space*/
											empty_pos[next_index*N*N + empty_index] = r*N + c;
											empty_index++;
										}
									}
								}
								new_array[next_index*N*N + row*N + col] = num;
								/*Record the number of empty spaces*/
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
	
	/*maximum number of boards from BFS*/
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
	
	for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            fscanf(fd, "%d", &Input[i * N + j]);
            //printf("%d\n", Input[i * N + j]);
            //if (!fscanf(a_file, "%c\n", &temp)) {
                //printf("File loading error!\n");
                //return;
            //}

            //if (temp >= '1' && temp <= '9') {
                //board[i * N + j] = (int) (temp - '0');
            //} else {
                //board[i * N + j] = 0;
            //}
        }
    }
	
	/* Initialize */
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
	/*number of times to do BFS, and loop_times is an even number here*/
	Kernel1<<<NumBlock, Blocksize>>>(pre_array, new_array, num_array, counter, empty_pos, num_empty);
	cudaMemcpy(&num_array, counter, sizeof(int), cudaMemcpyDeviceToHost);
	/*Loop to do BFS and to find all emepty position in the borad. Then save all new boards*/	
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
	/*Backtracking to complete the board*/
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
            printf("--------------------------------------------------\n");
        }

        for (int m = 0; m < N; m++) {
            if (m % n == 0) {
               printf("| ");
            }
            printf("%d ", Output[i * N + m]);
        }

        printf("|\n");
    }
    printf("--------------------------------------------------\n");
	printf("\n\nNow start to check the result\n");
    printf("------------------------------------------------------\n\n");
    /*check the result*/
	check(Output);
    /*free memory*/
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