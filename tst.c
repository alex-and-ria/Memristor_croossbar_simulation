




































#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include<string.h>
#include "omp.h"
#include"rw_.h"
//#include"node_schr.h"
#include"star_mesh_base_.c"

#include <unistd.h> // Required for sleep()


int main(int argc, char**argv){
     unsigned int** rw_; unsigned int** cl_; double** vl_;unsigned int* ln_;
     enum debug {mode1,mode2}; enum debug dbg;
     if(argc==3 && strcmp(argv[1],"debug_mode1")==0){
     	dbg=mode1;
     	dense_rdct(row,&rw_,col,&cl_,val,&vl_, &len1,&ln_, nds_td, &nds_n,1.,nds_td1,nds_n1, dbg,atoi(argv[2]));
     
     }
     else if(argc==3 && strcmp(argv[1],"debug_mode2")==0){
     	dbg=mode2;
     	dense_rdct(row,&rw_,col,&cl_,val,&vl_, &len1,&ln_, nds_td, &nds_n,1.,nds_td1,nds_n1, dbg,atoi(argv[2]));
     
     }
     else if(argc==2){//mode3 output;
     	dense_rdct(row,&rw_,col,&cl_,val,&vl_, &len1,&ln_, nds_td, &nds_n,atof(argv[1]),nds_td1,nds_n1,-1,0);
     
     }
     else{
     	printf("\nno mode found");
     	return 1;
     
     }
          //printf("%lld,%lld,%lld\n",mode1,mode2,mode1+mode2);
     

     return 0;

}
