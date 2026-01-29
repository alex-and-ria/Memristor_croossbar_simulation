




































#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "omp.h"
#include"rw_.h"
//#include"node_schr.h"
#include"star_mesh_base_.c"

#include <unistd.h> // Required for sleep()


int main(int argc, char**argv){
     unsigned int** rw_; unsigned int** cl_; double** vl_;unsigned int* ln_;
     //unsigned int m_n=4;
     //long long int mode1; long long int mode2;
     
          dense_rdct(row,&rw_,col,&cl_,val,&vl_, &len1,&ln_, nds_td, &nds_n,atof(argv[1]),nds_td1,nds_n1);
          //printf("%lld,%lld,%lld\n",mode1,mode2,mode1+mode2);
     

     return 0;

}
