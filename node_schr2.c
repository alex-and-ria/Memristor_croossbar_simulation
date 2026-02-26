
































//#define th_nb_koef (0.8)
#define fn "/share/share0/abc/star_mesh/log"

#include <time.h>
//#include<inttypes.h>
/*struct timespec curr_time;
ret=clock_gettime(CLOCK_MONOTONIC,&curr_time);
uint64_t tick = curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec;
mode1+=tick1-tick;*/


void node_analyzer(unsigned int* row, unsigned int *col, unsigned int *len, unsigned int* nds_td, unsigned int* nds_n, double* nb_koef, unsigned int* n_indp){//adj_m starts count from node 1 (Matlab indexing); it is represented in sparse Matlab format via [row,col,val]; col are sorted; nds_td should count from one (Matlab idexing), and be sorted; it is a list of nodes to delete; check nodes and delete one of them from nds_td if two are neighbours;
//TODO set data to file;
     unsigned int cnt_zr=0;
     unsigned int j=0;
     for(unsigned int i=0;i<*nds_n;i++){
          if(nds_td[i]==0) continue;
          while(col[j]<nds_td[i]) j++;
          
          for(unsigned int k=i+1; j<*len && k<*nds_n && col[j]==nds_td[i];){
               if(row[j]<nds_td[k]){
                    j++;
               }
               else if(row[j]>nds_td[k]){
                    k++;
               }
               else{
                    nds_td[k]=0;
                    j++; k++;
                    cnt_zr++;
               }
          }
          
     }
     (*nb_koef)=(cnt_zr+0.)/(*nds_n+0.);
     *n_indp=(*nds_n-cnt_zr);
}


void mode1_f(unsigned int **rw0,unsigned int** rw, unsigned int **cl0,unsigned int** cl, double **vl0, double** vl, unsigned int *len0, unsigned int *nds_td, unsigned int *nds_n,unsigned int **nds_td0,double th_nb_koef,unsigned char fl_dbg, unsigned int num_iter){
     double nb_koef=0; unsigned int n_indp=0;
     for(unsigned int i=0;i<*nds_n;i++){
          (*nds_td0)[i]=nds_td[i];
     }
     node_analyzer(*rw0,*cl0,len0,nds_td,nds_n, &nb_koef,&n_indp);
     
     
	unsigned int nds_n_nz=0; unsigned int nds_td_cnt=0; unsigned int *tmp_p_ui; double *tmp_p_d; unsigned int ln; unsigned int iter_cnt=0;
     while ((fl_dbg==0 && nb_koef<th_nb_koef && n_indp>1) || (fl_dbg==1 && iter_cnt<num_iter && n_indp>1)){//continue mode 1 until threshold or until number of indipendent nodes more then 1;
          nds_n_nz=0; nds_td_cnt=0;
          for(unsigned int i=0;i<*nds_n;){//delete zeros; keep nodes to delete on next iterations;
               while(i<*nds_n && nds_td[i]==0){
                    nds_td[nds_td_cnt]=(*nds_td0)[i];
                    nds_td_cnt++;
                    i++;
               
               }
               while(i<*nds_n && nds_td[i]!=0){
                    (*nds_td0)[nds_n_nz]=(*nds_td0)[i];
                    nds_n_nz++;
                    i++;
                    
               }
               
          }
          //at this point nds_td should contain list of nodes that are not neighbours, and, hece, safe to delete;
          *nds_n=nds_n_nz;
          star_mesh_base(*rw0,rw,*cl0,cl,*vl0,vl,len0,&ln,(*nds_td0),nds_n);
          
          
          tmp_p_ui=(*cl); (*cl)=(*cl0); (*cl0)=tmp_p_ui; free(*cl);
          tmp_p_ui=(*rw); (*rw)=(*rw0); (*rw0)=tmp_p_ui; free(*rw);
          tmp_p_d=(*vl); (*vl)=(*vl0); (*vl0)=tmp_p_d; free(*vl);
          
          
          (*len0)=ln;
          *nds_n=nds_td_cnt;//since number of nodes to delete is decreesing at each new iteration, no need to realloc memory;
          for(unsigned int i=0;i<*nds_n;i++){
               (*nds_td0)[i]=nds_td[i];
          }
          node_analyzer(*rw0,*cl0,len0,nds_td,nds_n, &nb_koef,&n_indp);
          iter_cnt++;
          
     
     
     }

}



void mode2_f(unsigned int **rw0,unsigned int **rw, unsigned int **cl0, unsigned int **cl, double **vl0, double **vl, unsigned int *len0, unsigned int *nds_td0, unsigned int nds_n,unsigned char fl_dbg, unsigned int num_iter,unsigned char fl_parallel){
     unsigned int *tmp_p_ui; double *tmp_p_d;
     unsigned int* rw_msh=(unsigned int*)malloc(((*len0)*((*len0)-1))*sizeof(unsigned int));//allocate space for maximum number of new edges (if all nodes were connected); potentially can reallocate each time for (j1-j0)*(j1-j0-1) each iteration for memory usage optimization;
     unsigned int* cl_msh=(unsigned int*)malloc(((*len0)*((*len0)-1))*sizeof(unsigned int));
     double* vl_msh=(double*)malloc(((*len0)*((*len0)-1))*sizeof(double));
     double cnds_sum=0; unsigned int cnt_curr=0;
     unsigned int cnt_m=0; unsigned int iter_cnt=0;
     
     for(unsigned int i=0,j=0,j0=0,j1=0;i<*len0 && j<nds_n && (fl_dbg==0 || (fl_dbg==1 && iter_cnt<num_iter));){//*nds_n is smaller or equal to *len0;
          while((*cl0)[i]<nds_td0[j]){
               i++;
          
          }
          j0=i;
          cnds_sum=0;
          while(i<*len0 && (*cl0)[i]==nds_td0[j]){
               cnds_sum+=(*vl0)[i];
               i++;
               
          }
          j1=i;//index boundearies for current node to delete;
          
          ////////////////
          #ifdef _OPENMP
               unsigned int num_threads_val=((j1-1-j0)<(unsigned int)omp_get_max_active_levels())?(j1-1-j0):(unsigned int) omp_get_max_active_levels();
          #endif
          #pragma omp parallel for num_threads(num_threads_val) lastprivate(cnt_curr) if(fl_parallel==1)
          for(unsigned int kk=j0;kk<j1-1;kk++){
               cnt_curr=0;
               
               for(unsigned int ll=kk+1;ll<j1;ll++){//TODO_ check;
                    cnt_curr=(j1-j0-1)*(kk-j0)+(ll-j0-1);
                    cl_msh[cnt_curr]=(*rw0)[kk];//(*cl) gets smaller index;
                    rw_msh[cnt_curr]=(*rw0)[ll];//(*rw) gets bigger index;
                    vl_msh[cnt_curr]=1./(((1./(*vl0)[kk])*(1./(*vl0)[ll]))*cnds_sum);
                    
                    cnt_curr=(j1-j0-1)*(ll-j0)+(kk-j0);
                    cl_msh[cnt_curr]=(*rw0)[ll];//add duplicated to comply with adjacency matrix format;
                    rw_msh[cnt_curr]=(*rw0)[kk];
                    vl_msh[cnt_curr]=1./(((1./(*vl0)[kk])*(1./(*vl0)[ll]))*cnds_sum);
                    
               
               }
               
          
          }
          //////////
          
          cnt_curr++;
          unsigned int k0=0, k1=0; cnt_m=0;
          (*rw)=(unsigned int*)malloc(((j1-j0)*(j1-j0-1)+(*len0))*sizeof(unsigned int));//allocate space for number of new edges in the resulting (after node delteion) mesh and all the remaining nodes;
          (*cl)=(unsigned int*)malloc(((j1-j0)*(j1-j0-1)+(*len0))*sizeof(unsigned int));
          (*vl)=(double*)malloc(((j1-j0)*(j1-j0-1)+(*len0))*sizeof(double));
          for(;k0<*len0 && k1<cnt_curr;){
          	if(k0==j0) k0=j1;
          	if((*rw0)[k0]==nds_td0[j]){k0++; continue;}
          	if(((*cl0)[k0]<cl_msh[k1]) || ((*cl0)[k0]==cl_msh[k1] && (*rw0)[k0]<rw_msh[k1])){
          		(*cl)[cnt_m]=(*cl0)[k0];
          		(*rw)[cnt_m]=(*rw0)[k0];
          		(*vl)[cnt_m]=(*vl0)[k0];
          		cnt_m++; k0++;
          		
          	
          	}
          	else if((*cl0)[k0]==cl_msh[k1] && (*rw0)[k0]==rw_msh[k1]){
          		(*cl)[cnt_m]=(*cl0)[k0];
          		(*rw)[cnt_m]=(*rw0)[k0];
          		(*vl)[cnt_m]=(*vl0)[k0]+vl_msh[k1];
          		cnt_m++; k0++; k1++; 
          	
          	}
          	else{
          		(*cl)[cnt_m]=cl_msh[k1];
          		(*rw)[cnt_m]=rw_msh[k1];
          		(*vl)[cnt_m]=vl_msh[k1];
          		cnt_m++; k1++;
          	
          	}
          	
          
          }
          for(;k0<*len0;){
          	if(k0==j0) k0=j1;
          	if((*rw0)[k0]==nds_td0[j]){k0++; continue;}
     		(*cl)[cnt_m]=(*cl0)[k0];
     		(*rw)[cnt_m]=(*rw0)[k0];
     		(*vl)[cnt_m]=(*vl0)[k0];
     		cnt_m++; k0++;
          
          
          }
          for(;k1<cnt_curr;){
			(*cl)[cnt_m]=cl_msh[k1];
			(*rw)[cnt_m]=rw_msh[k1];
			(*vl)[cnt_m]=vl_msh[k1];
			cnt_m++; k1++;
          		
          
          }
          
          j++; i=0; iter_cnt++;
          
          tmp_p_ui=(*cl0); (*cl0)=(*cl); free(tmp_p_ui);
          tmp_p_ui=(*rw0); (*rw0)=(*rw); free(tmp_p_ui);
          tmp_p_d=(*vl0); (*vl0)=(*vl); free(tmp_p_d);
          (*rw0)=(unsigned int*)realloc((*rw0),cnt_m*sizeof(unsigned int));
          (*cl0)=(unsigned int*)realloc((*cl0),cnt_m*sizeof(unsigned int));
          (*vl0)=(double*)realloc((*vl0),cnt_m*sizeof(double));
          *len0=cnt_m;
     
     
     }
     free(rw_msh); free(cl_msh); free(vl_msh);
     

}

void mode3_f(unsigned int ***rw0,unsigned int *rw00, unsigned int ***cl0, unsigned int *cl00, double ***vl0, double *vl00, unsigned int **len0,unsigned int ln, unsigned int *nds_td1, unsigned int nds_n1, unsigned int max_m_sz, unsigned int* n_th){
     (*n_th)=(nds_n1%max_m_sz==0)?nds_n1/max_m_sz:nds_n1/max_m_sz+1;
     unsigned int **nds_td0=(unsigned int**)malloc((*n_th)*sizeof(unsigned int*));
     (*cl0)=(unsigned int**)malloc((*n_th)*sizeof(unsigned int*));
     (*rw0)=(unsigned int**)malloc((*n_th)*sizeof(unsigned int*));
     (*vl0)=(double**)malloc((*n_th)*sizeof(double*));
     unsigned int **cl=(unsigned int**)malloc((*n_th)*sizeof(unsigned int*));
     unsigned int **rw=(unsigned int**)malloc((*n_th)*sizeof(unsigned int*));
     double **vl=(double**)malloc((*n_th)*sizeof(double*));
     (*len0)=(unsigned int*)malloc((*n_th)*sizeof(unsigned int));

     #ifdef _OPENMP
          unsigned int num_threads_val=((*n_th)<(unsigned int)omp_get_max_active_levels())?(*n_th):(unsigned int)omp_get_max_active_levels();
     #endif
     #pragma omp parallel for num_threads(num_threads_val)
     //parallel;
     for(unsigned int i=0;i<(*n_th);i++){
          (*len0)[i]=ln;
          (*cl0)[i]=(unsigned int*)malloc((*len0)[i]*sizeof(unsigned int));
          (*rw0)[i]=(unsigned int*)malloc((*len0)[i]*sizeof(unsigned int));
          (*vl0)[i]=(double*)malloc((*len0)[i]*sizeof(double));
          nds_td0[i]=(unsigned int*)malloc((nds_n1-max_m_sz)*sizeof(unsigned int));//all threads except the last one will keep max_m_sz nodes (so delete other nodes);
          for(unsigned int j=0;j<(*len0)[i];j++){
               (*cl0)[i][j]=cl00[j];
               (*rw0)[i][j]=rw00[j];
               (*vl0)[i][j]=vl00[j];
          
          }
     /*}NOTE: merged to loops together for openmp efficiency;
     if(nds_n1%max_m_sz!=0) {nds_td0[(*n_th)-1]=(unsigned int*)realloc(nds_td0[(*n_th)-1],(nds_n1-nds_n1%max_m_sz)*sizeof(unsigned int));}//this thread should delete all nodes except ther remaining;
     
     //parallel;
     for(unsigned int i=0;i<(*n_th);i++){*/
          unsigned int l_idx=i*max_m_sz, r_idx=((i+1)*max_m_sz<=nds_n1)?(i+1)*max_m_sz:nds_n1;
          unsigned int j=0;
          for(;j<l_idx;j++){
               nds_td0[i][j]=nds_td1[j];
          
          }
          for(j=r_idx;j<nds_n1;j++){
               nds_td0[i][j-(r_idx-l_idx)]=nds_td1[j];
          
          }
          //after mode1 the number of indipendent nodes reduced, after mode2 graph become even denser, so at this point mostly one by one node deletion shoud be performed;
          mode2_f(&((*rw0)[i]),&(rw[i]),&((*cl0)[i]),&(cl[i]),&((*vl0)[i]),&(vl[i]), &((*len0)[i]),nds_td0[i],(nds_n1-(r_idx-l_idx)),0,0,0);
          //mode2_f(unsigned int **rw0,unsigned int **rw, unsigned int **cl0, unsigned int **cl, double **vl0, double **vl, unsigned int *len0, unsigned int *nds_td0, unsigned int nds_n,unsigned char fl_dbg, unsigned int num_iter);
          free(nds_td0[i]);
     
     }
     
     free(cl); free(rw); free(vl);
     free(nds_td0);
     

}

void dense_rdct(unsigned int *row, unsigned long long int* rw_v, unsigned int *col, unsigned long long int* cl_v, double *val, unsigned long long int* vl_v, unsigned int *len,unsigned int **ln_, unsigned int *nds_td, unsigned int *nds_n,double th_nb_koef, unsigned int *nds_td1, unsigned int nds_n1, unsigned int* n_th,unsigned int max_m_sz, int mode_dbg, unsigned int num_iter,unsigned int m_n){//unsigned long long int* here acts as generic void*, but stored as plain 64-bit number;
	enum debug {mode1,mode2};
	unsigned int*** rw_=(unsigned int***) rw_v; unsigned int*** cl_=(unsigned int***) cl_v; double*** vl_=(double***) vl_v;
     unsigned int *nds_td0=(unsigned int*)malloc((*nds_n)*sizeof(unsigned int)); //with more nodes deleted, graph becomes more connected, so less indipendent (that are not neighboues), hence initial memory allocation here should suffice;
     unsigned int len0=*len;
     
     unsigned int *rw0=(unsigned int*)malloc(len0*sizeof(unsigned int));
     unsigned int *cl0=(unsigned int*)malloc(len0*sizeof(unsigned int));
     double *vl0=(double*)malloc(len0*sizeof(double));
     for(unsigned int i=0;i<len0;i++){
          cl0[i]=col[i]; rw0[i]=row[i]; vl0[i]=val[i];
     
     }
     char fn_buff[128]; sprintf(fn_buff,"%dx%d.csv",m_n,m_n);
     FILE* fp=fopen (fn_buff,"ab");
     struct timespec curr_time;
     long long unsigned int tick;
     long long unsigned int tick1;
     long long unsigned int mode1_time=0, mode2_time=0,mode3_time=0,tot_time=0;
     fprintf(fp,"\nth_nb_koef:,%f",th_nb_koef);
     
     
     if(mode_dbg==-1 && max_m_sz<nds_n1){
          unsigned int* rw=NULL; unsigned int* cl=NULL;double* vl=NULL;
     	clock_gettime(CLOCK_MONOTONIC,&curr_time); tick = curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec;
     	mode1_f(&rw0,&rw,&cl0,&cl,&vl0,&vl,&len0,nds_td,nds_n,&nds_td0,th_nb_koef,0,0);//input is rw0,cl0,vl0 of size len0; output is rw0,cl0,vl0, of size len0, (inout parameters), nds_td and nds_td0 are of size *nds_n (inout), th_nb_koef is threshold parameter (input), mode1 is time spent in mode 1 (ouput);
     	clock_gettime(CLOCK_MONOTONIC,&curr_time); tick1 = curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec;
     	mode1_time=tick1-tick;
     	
     	clock_gettime(CLOCK_MONOTONIC,&curr_time); tick = curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec;
     	mode2_f(&rw0,&rw,&cl0,&cl,&vl0,&vl,&len0,nds_td0,*nds_n,0,0,1);//input is rw0,cl0,vl0 of size len0; output is rw0,cl0,vl0, of size len0, (inout parameters), nds_td0 is of size *nds_n (input), mode2 is time spent in mode2 (output);
     	clock_gettime(CLOCK_MONOTONIC,&curr_time); tick1 = curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec;
     	mode2_time=tick1-tick;
     	
     	clock_gettime(CLOCK_MONOTONIC,&curr_time); tick = curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec;
     	mode3_f(rw_,rw0,cl_,cl0, vl_,vl0,ln_,len0, nds_td1,nds_n1, max_m_sz,n_th);//input is in rw0,cl0,vl0,len0,nds_td1, nds_n1,max_m_sz, inputs are not modified here; the rest of the parameters are outputs; n_th is a number of threads; rw_,cl_, and vl_ should reutrn *n_th arrays, each containing correspoinding rows,cols, and vals of corresponding matrixes; ln_ should return n_th sizes of resulting rw_,cl_, and vl_ arrays;
     	//mode3_f(unsigned int ***rw0,unsigned int *rw00, unsigned int ***cl0, unsigned int *cl00, double ***vl0, double *vl00, unsigned int **len0,unsigned int ln, unsigned int *nds_td1, unsigned int nds_n1, unsigned int max_m_sz, unsigned int* n_th)
     	clock_gettime(CLOCK_MONOTONIC,&curr_time); tick1 = curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec;
     	mode3_time=tick1-tick;
     	
     	free(rw0); free(cl0); free(vl0);
     
     }
     else if(mode_dbg==-1 && max_m_sz>=nds_n1){
          (*rw_)=(unsigned int**)malloc(1*sizeof(unsigned int*));// (*rw_)[0]=(unsigned int*)malloc(1*sizeof(unsigned int*));
		(*cl_)=(unsigned int**)malloc(1*sizeof(unsigned int*));// (*cl_)[0]=(unsigned int*)malloc(1*sizeof(unsigned int*));
		(*vl_)=(double**)malloc(1*sizeof(double*));// (*vl_)[0]=(double*)malloc(1*sizeof(double*));
		(*ln_)=(unsigned int*)malloc(1*sizeof(unsigned int));
		clock_gettime(CLOCK_MONOTONIC,&curr_time); tick = curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec;
		mode1_f(&rw0,&((*rw_)[0]),&cl0,&((*cl_)[0]),&vl0,&((*vl_)[0]),&len0,nds_td,nds_n,&nds_td0,th_nb_koef,0,0);
		clock_gettime(CLOCK_MONOTONIC,&curr_time); tick1 = curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec;
     	mode1_time=tick1-tick;
     	
     	clock_gettime(CLOCK_MONOTONIC,&curr_time); tick = curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec;
		mode2_f(&rw0,&((*rw_)[0]),&cl0,&((*cl_)[0]),&vl0,&((*vl_)[0]),&len0,nds_td0,*nds_n,0,0,1);
		clock_gettime(CLOCK_MONOTONIC,&curr_time); tick1 = curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec;
     	mode2_time=tick1-tick;
		(*rw_)[0]=rw0; (*cl_)[0]=cl0; (*vl_)[0]=vl0;
		(*ln_)[0]=len0;
		(*n_th)=1;
     
     
     }
     else if(mode_dbg==mode1){
     	(*rw_)=(unsigned int**)malloc(1*sizeof(unsigned int*)); //(*rw_)[0]=(unsigned int*)malloc(1*sizeof(unsigned int*));
		(*cl_)=(unsigned int**)malloc(1*sizeof(unsigned int*)); //(*cl_)[0]=(unsigned int*)malloc(1*sizeof(unsigned int*));
		(*vl_)=(double**)malloc(1*sizeof(double*)); //(*vl_)[0]=(double*)malloc(1*sizeof(double*));
		(*ln_)=(unsigned int*)malloc(1*sizeof(unsigned int));
     	mode1_f(&rw0,&((*rw_)[0]),&cl0,&((*cl_)[0]),&vl0,&((*vl_)[0]),&len0,nds_td,nds_n,&nds_td0,th_nb_koef,1,num_iter);
     	(*rw_)[0]=rw0; (*cl_)[0]=cl0; (*vl_)[0]=vl0;
     	(*ln_)[0]=len0;
     	(*n_th)=1;
     
     }
     else if(mode_dbg==mode2){
     	(*rw_)=(unsigned int**)malloc(1*sizeof(unsigned int*));// (*rw_)[0]=(unsigned int*)malloc(1*sizeof(unsigned int*));
		(*cl_)=(unsigned int**)malloc(1*sizeof(unsigned int*));// (*cl_)[0]=(unsigned int*)malloc(1*sizeof(unsigned int*));
		(*vl_)=(double**)malloc(1*sizeof(double*));// (*vl_)[0]=(double*)malloc(1*sizeof(double*));
		(*ln_)=(unsigned int*)malloc(1*sizeof(unsigned int));
		mode1_f(&rw0,&((*rw_)[0]),&cl0,&((*cl_)[0]),&vl0,&((*vl_)[0]),&len0,nds_td,nds_n,&nds_td0,th_nb_koef,0,0);
		mode2_f(&rw0,&((*rw_)[0]),&cl0,&((*cl_)[0]),&vl0,&((*vl_)[0]),&len0,nds_td0,*nds_n,1,num_iter,1);
		(*rw_)[0]=rw0; (*cl_)[0]=cl0; (*vl_)[0]=vl0;
		(*ln_)[0]=len0;
		(*n_th)=1;
     
     }
     tot_time=mode1_time+mode2_time+mode3_time;
     fprintf(fp,",%llu",mode1_time);
     fprintf(fp,",%llu",mode2_time);
     fprintf(fp,",%llu",mode3_time);
     fprintf(fp,",%llu",tot_time);
     fclose(fp);
     
     
     //(*cl)=cl0; (*rw)=rw0; (*vl)=vl0;
     //*ln=len0;
    
     
     free(nds_td0);

 
     //fclose(fp);


}

void get_dbg_arr(unsigned long long int *rw_v,unsigned long long int* cl_v, unsigned long long int* vl_v, unsigned int** rw0, unsigned int** cl0, double **vl0, unsigned int iter){
     (*rw0)=((unsigned int**) (*rw_v))[iter];
     (*cl0)=((unsigned int**) (*cl_v))[iter];
     (*vl0)=((double**) (*vl_v))[iter];

}


void data_free(unsigned long long int *rw_v,unsigned long long int* cl_v, unsigned long long int* vl_v, unsigned int **ln,unsigned int *n_th){
     unsigned int** rw=(unsigned int**) (*rw_v);
     unsigned int** cl=(unsigned int**) (*cl_v);
     double** vl=(double**) (*vl_v);
     for(unsigned int i=0;i<(*n_th);i++){
          free(rw[i]); free(cl[i]); free(vl[i]);
     
     }
     free(rw); free(cl); free(vl);
     free(*ln);

}

/*void get_pointer(long long unsigned int* ptr_tmp, unsigned int num_out){
     (*ptr_tmp)= (long long unsigned int) malloc(num_out*sizeof(unsigned int*));

}*/

