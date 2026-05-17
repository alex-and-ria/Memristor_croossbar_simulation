




































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
     enum debug {mode1,mode2,mode_1_2_3}; enum debug dbg;
     unsigned int n_th;//, max_m_sz=50;
     if(argc==3 && strcmp(argv[1],"debug_mode1")==0){
     	dbg=mode1;
     	dense_rdct(row,(unsigned long long int*)(&rw_),col,(unsigned long long int*)(&cl_),val,(unsigned long long int*)(&vl_), &len1,&ln_, nds_td, &nds_n,1.,nds_td1,nds_n1, &n_th,max_m_sz,dbg,atoi(argv[2]));
     	//dense_rdct(unsigned int *row, unsigned long long int* rw_v, unsigned int *col, unsigned long long int* cl_v, double *val, unsigned long long int* vl_v, unsigned int *len,unsigned int **ln_, 
     	     //unsigned int *nds_td, unsigned int *nds_n,double th_nb_koef, unsigned int *nds_td1, unsigned int nds_n1, unsigned int* n_th,unsigned int max_m_sz, int mode_dbg, unsigned int num_iter)
     
     }
     else if(argc==3 && strcmp(argv[1],"debug_mode2")==0){
     	dbg=mode2;
     	dense_rdct(row,(unsigned long long int*)(&rw_),col,(unsigned long long int*)(&cl_),val,(unsigned long long int*)(&vl_), &len1,&ln_, nds_td, &nds_n,1.,nds_td1,nds_n1, &n_th,max_m_sz, dbg,atoi(argv[2]));
     
     }
     else if(argc==2 &&0){//mode3 output;
     	dense_rdct(row,(unsigned long long int*)(&rw_),col,(unsigned long long int*)(&cl_),val,(unsigned long long int*)(&vl_), &len1,&ln_, nds_td, &nds_n,atof(argv[1]),nds_td1,nds_n1,&n_th,max_m_sz,-1,0);
     
     }
     else if(argc==1 || argc==2){
          //star_mesh_base_alg(row,rw_,col,cl_,val,vl_, len1,ln_, nds_td, &nds_n);
          //star_mesh_base_alg(unsigned int *row,unsigned int** rw, unsigned int *col,unsigned int** cl, double *val,/*double** vl,*/ unsigned int len/*,unsigned int *ln, unsigned int *nds_td, unsigned int *nds_n*/){
          /*unsigned char* fl_nm; unsigned int** tot_mem_buff; unsigned int mem_buff_cnt;
          schr_rec(row,col,len1,nds_td,nds_n,0.75,4,4,&fl_nm,&tot_mem_buff,&mem_buff_cnt);
          //schr_rec(unsigned int *rw,unsigned int *cl, unsigned int len, unsigned int *nds_td, unsigned int nds_n,double th_nb_koef,unsigned int m, unsigned int n, unsigned char** fl_nm, unsigned int*** tot_mem_buff, unsigned int* mem_buff_cnt)
          schr_rec_free(tot_mem_buff,mem_buff_cnt,fl_nm);
          //schr_rec_free(unsigned int*** tot_mem_buff, unsigned int* mem_buff_cnt,unsigned char** fl_nm)
          schr_read(nds_td, nds_n,0.75,4,4,row,col,val,len1);
          //schr_read(unsigned int *nds_td, unsigned int nds_n,double th_nb_koef,unsigned int m, unsigned int n, unsigned int *row, unsigned int *col, unsigned int* val,unsigned int len)*/
          
     }
     else if(argc==4 && strcmp(argv[3],"debug_mode_alg")==0){
          dbg=mode_1_2_3;
          dense_rdct(row,(unsigned long long int*)(&rw_), col, (unsigned long long int*)(&cl_), val, (unsigned long long int*)(&vl_), &len1,&ln_, nds_td, &nds_n,atof(argv[1]), nds_td1, nds_n1, &n_th,max_m_sz, mode_1_2_3,0);
     
     }
     else{
     	printf("\nno mode found");
     	return 1;
     
     }
          //printf("%lld,%lld,%lld\n",mode1,mode2,mode1+mode2);
     

     return 0;

}
