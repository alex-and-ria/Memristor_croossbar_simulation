#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>
int main(){
printf("sz=%ld\nlog of size_max: %lf\nsizeof(size_t)=%ld\n",sizeof(unsigned long long int),log(SIZE_MAX)/log(2),sizeof(size_t)); return 0;


}
