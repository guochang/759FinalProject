#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define n 3

void main(){
	int a = n*(8/n);
	int b = 8 / n;
	int c = n*b;
	printf("%d,%d,%d",a,b,c);
}