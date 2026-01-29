
































//#define th_nb_koef (0.8)
#define fn "/share/share0/abc/star_mesh/log"

#include <time.h>
//#include<inttypes.h>
/*struct timespec curr_time;
ret=clock_gettime(CLOCK_MONOTONIC,&curr_time);
uint64_t tick = curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec;*/


void node_analyzer(unsigned int* row, unsigned int *col, unsigned int *len, unsigned int* nds_td, unsigned int* nds_n, double* nb_koef, unsigned int* n_indp){//adj_m starts count from node 1 (Matlab indexing); it is represented in sparse Matlab format via [row,col,val]; col are sorted; nds_td should count from one (Matlab idexing), and be sorted; it is a list of nodes to delete; check nodes and delete one of them from nds_td if two are neighbours;
//TODO set data to file;
     //unsigned int *idx_l, *idx_r, idxl,idxr;  idx_l=&idxl; idx_r=&idxr;
     unsigned int cnt_zr=0;
     unsigned int j=0;
     for(unsigned int i=0;i<*nds_n;i++){
          if(nds_td[i]==0) continue;
          //col_fnd(col,len, nds_td[i],idx_l,idx_r);
          while(col[j]<nds_td[i]) j++;
          
          for(unsigned int /*j=*idx_l,*/k=i+1; j<*len && k<*nds_n && col[j]==nds_td[i];){
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

void dense_rdct(unsigned int *row,unsigned int** rw, unsigned int *col,unsigned int** cl, double *val, double** vl, unsigned int *len,unsigned int *ln, unsigned int *nds_td, unsigned int *nds_n,double th_nb_koef,long long int* mode1, long long int* mode2){
     //char fn_buff[128]; sprintf(fn_buff,"%s_%dx%d_%f.csv",fn,m_n,m_n,th_nb_koef);
     //FILE* fp=fopen (fn_buff,"wb");
     struct timespec curr_time;
     long long unsigned int tick;
     long long unsigned int tick1;
     *mode1=0; *mode2=0;
     
     double nb_koef=0; unsigned int n_indp=0;
     unsigned int *nds_td0=(unsigned int*)malloc((*nds_n)*sizeof(unsigned int)); unsigned int nds_td_cnt=0;//with more nodes deleted, graph becomes more connected, so less indipendent (that are not neighboues), hence initial memory allocation here should suffice;
     //unsigned int *nds_td1=(unsigned int*)malloc((*nds_n)*sizeof(unsigned int));//initial memory allocation here should suffice;
     for(unsigned int i=0;i<*nds_n;i++){
          nds_td0[i]=nds_td[i];
     }
     clock_gettime(CLOCK_MONOTONIC,&curr_time); tick = curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec;
     node_analyzer(row,col,len,nds_td,nds_n, &nb_koef,&n_indp);
     clock_gettime(CLOCK_MONOTONIC,&curr_time); tick1 = curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec;
     //fprintf(fp,"\nn_indp:,%u",n_indp);
     //fprintf(fp,",node_analyzer:,%llu",tick1-tick);
     *mode1+=tick1-tick;
     
     
     unsigned int nds_n_nz=0;
     //unsigned int swp_tmp=0; double swp_tmp_d=0;
     /*(*rw)=(unsigned int*)malloc((*len)*sizeof(unsigned int));
     (*cl)=(unsigned int*)malloc((*len)*sizeof(unsigned int));
     (*vl)=(double*)malloc((*len)*sizeof(double));*/
     unsigned int *rw0=(unsigned int*)malloc((*len)*sizeof(unsigned int));
     unsigned int *cl0=(unsigned int*)malloc((*len)*sizeof(unsigned int));
     double *vl0=(double*)malloc((*len)*sizeof(double));
     unsigned int *tmp_p_ui; double *tmp_p_d;
     unsigned int len0=*len;
     for(unsigned int i=0;i<*len;i++){
          cl0[i]=col[i]; rw0[i]=row[i]; vl0[i]=val[i];
     
     }
     //rw0=row; cl0=col; vl0=val;
     while (nb_koef<th_nb_koef && n_indp>1){//continue mode 1 until threshold or until number of indipendent nodes more then 1;
          nds_n_nz=0; nds_td_cnt=0;
          for(unsigned int i=0;i<*nds_n;){//delete zeros; keep nodes to delete on next iterations;
               while(i<*nds_n && nds_td[i]==0){
                    nds_td[nds_td_cnt]=nds_td0[i];
                    nds_td_cnt++;
                    i++;
               
               }
               while(i<*nds_n && nds_td[i]!=0){
                    nds_td0[nds_n_nz]=nds_td0[i];
                    nds_n_nz++;
                    i++;
                    
               }
               
          }
          //at this point nds_td should contain list of nodes that are not neighbours, and, hece, safe to delete;
          *nds_n=nds_n_nz;
          /*if(*nds_n==0) break;
          if(*nds_n==1) {
               *nds_n=nds_td_cnt;//since number of nodes to delete is decreesing at each new iteration, no need to realloc memory;
               for(unsigned int i=0;i<*nds_n;i++){
                    nds_td0[i]=nds_td1[i]; nds_td[i]=nds_td1[i];
               }
               goto m;
          }*/
          //fprintf(fp,"\nnodes:,%d",*nds_n);
          clock_gettime(CLOCK_MONOTONIC,&curr_time); tick = curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec;
          star_mesh_base(rw0,rw,cl0,cl,vl0,vl,&len0,ln,nds_td0,nds_n);
          clock_gettime(CLOCK_MONOTONIC,&curr_time); tick1 = curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec;
          //fprintf(fp,",mode 1:,%llu",tick1-tick);
          *mode1+=tick1-tick;
          
          /*clock_gettime(CLOCK_MONOTONIC,&curr_time); tick = curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec;
          for(unsigned int i=0;i<*ln;i++){//sort with selection sort (major is column sorting, secondary is row sorting); faster alghourithms might be equivalently used (like mege sort, or quick sort);
               for(unsigned int j=i+1;j<*ln;j++){
                    if((*cl)[j]<(*cl)[i] || ((*cl)[j]==(*cl)[i] && (*rw)[j]<((*rw)[i]))){
                         swp_tmp=(*cl)[i]; (*cl)[i]=(*cl)[j]; (*cl)[j]=swp_tmp;
                         swp_tmp=(*rw)[i]; (*rw)[i]=(*rw)[j]; (*rw)[j]=swp_tmp;
                         swp_tmp_d=(*vl)[i]; (*vl)[i]=(*vl)[j]; (*vl)[j]=swp_tmp_d;
                     
                    }
               
               }
          
          }//at this point *rw,*cl, and *vl should be sorted like standart Matlab column major format;
          
          clock_gettime(CLOCK_MONOTONIC,&curr_time); tick1 = curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec;
          fprintf(fp,"\nsort:,%llu",tick1-tick);*/
          
          tmp_p_ui=(*cl); (*cl)=cl0; cl0=tmp_p_ui; free(*cl);
          tmp_p_ui=(*rw); (*rw)=rw0; rw0=tmp_p_ui; free(*rw);
          tmp_p_d=(*vl); (*vl)=vl0; vl0=tmp_p_d; free(*vl);
          
          
          len0=*ln;
          *nds_n=nds_td_cnt;//since number of nodes to delete is decreesing at each new iteration, no need to realloc memory;
          for(unsigned int i=0;i<*nds_n;i++){
               nds_td0[i]=nds_td[i];
          }
          clock_gettime(CLOCK_MONOTONIC,&curr_time); tick = curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec;
          node_analyzer(rw0,cl0,&len0,nds_td,nds_n, &nb_koef,&n_indp);
          clock_gettime(CLOCK_MONOTONIC,&curr_time); tick1 = curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec;
          //fprintf(fp,"\nn_indp:,%u",n_indp);
          //fprintf(fp,",node_analyzer:,%llu",tick1-tick);
          *mode1+=tick1-tick;
          
     
     
     
     }
     //m:
     //fprintf(fp,"\nmode 2:\n");
     //nds_td0 nds_n rw0 cl0 vl0 len
     //(*rw)=(unsigned int*)malloc(((*len)*((*len)-1))*sizeof(unsigned int));//allocate space for maximum number of new edges (if all nodes were connected);
     //(*cl)=(unsigned int*)malloc(((*len)*((*len)-1))*sizeof(unsigned int));
     //(*vl)=(double*)malloc(((*len)*((*len)-1))*sizeof(double));
     //unsigned int *nz_nd_a=(unsigned int*) malloc(2*((len0)-(*nds_n))*sizeof(unsigned int));//approximate number that is sufficient to keep reamaining nodes;
     unsigned int* rw_msh=(unsigned int*)malloc(((*len)*((*len)-1))*sizeof(unsigned int));//allocate space for maximum number of new edges (if all nodes were connected);
     unsigned int* cl_msh=(unsigned int*)malloc(((*len)*((*len)-1))*sizeof(unsigned int));
     double* vl_msh=(double*)malloc(((*len)*((*len)-1))*sizeof(double));
     double cnds_sum=0; unsigned int cnt_curr=0;//,cnt_curr0=0;//,cnt_curr_init=0; //unsigned int cnt_nz_nds=0; unsigned char fl_nz=0; //unsigned int nds_i_val=0;
     unsigned int cnt_m=0;
     for(unsigned int i=0,j=0,j0=0,j1=0;i<len0 && j<*nds_n;){//*nds_n is smaller or equal to len0;
          clock_gettime(CLOCK_MONOTONIC,&curr_time); tick = curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec;
          while(cl0[i]<nds_td0[j]){
               i++;
          
          }
          j0=i;
          cnds_sum=0;
          while(i<len0 && cl0[i]==nds_td0[j]){
               cnds_sum+=vl0[i];
               i++;
               
          }
          j1=i;//index boundearies for current node to delete;
          
          //unsigned int curr_sz=((j1-j0)*(j1-j0-1))/2.,cnt_curr0=0;
          ////////////////
          //parallel
          for(unsigned int kk=j0;kk<j1-1;kk++){
               cnt_curr=0;
               //unsigned int kk0_tmp=kk-j0;
               //cnt_curr=kk0_tmp*(j1-j0)-((1+kk0_tmp)*kk0_tmp)/2;//calculating index to avoid dependencies with previous iterations;  sum_{i=1}^{i=kk0_tmp} (j1-j0-i);
               /*for(unsigned int kk0=j1-j0-1;kk0!=j1-kk-1;kk0--){
                    cnt_curr+=kk0;
               
               }*///calculating index to avoid dependencies with previous iterations;
               
               for(unsigned int ll=kk+1;ll<j1;ll++){//TODO_ check;
                    cnt_curr=(j1-j0-1)*(kk-j0)+(ll-j0-1);
                    /*cl_msh[cnt_curr+ll-kk-1]=rw0[kk];//(*cl) gets smaller index;
                    rw_msh[cnt_curr+ll-kk-1]=rw0[ll];//(*rw) gets bigger index;
                    vl_msh[cnt_curr+ll-kk-1]=1./((1./vl0[kk])*(1./vl0[ll])*cnds_sum);*/
                    cl_msh[cnt_curr]=rw0[kk];//(*cl) gets smaller index;
                    rw_msh[cnt_curr]=rw0[ll];//(*rw) gets bigger index;
                    vl_msh[cnt_curr]=1./((1./vl0[kk])*(1./vl0[ll])*cnds_sum);
                    //cnt_curr0=(j1-j0-1)*(ll-j0)+(kk-j0);
                    cnt_curr=(j1-j0-1)*(ll-j0)+(kk-j0);
                    //cnt_curr0=curr_sz+((1+(ll-j0-1))*(ll-j0-1))/2 + (kk-j0);//for second index, we start after curr_sz (size of first part); we count how many elements are from 1 to (ll-j0-1), then adjust index by (kk-j0);
                    cl_msh[cnt_curr]=rw0[ll];//add duplicated to comply with adjacency matrix format;
                    rw_msh[cnt_curr]=rw0[kk];
                    vl_msh[cnt_curr]=vl_msh[cnt_curr];
                    
               
               }
               
          
          }
          //////////
          //cnt_curr_init=cnt_curr+1;
          cnt_curr++;//cnt_curr=2*curr_sz;
          unsigned int k0=0, k1=0; cnt_m=0;
          (*rw)=(unsigned int*)malloc(((j1-j0)*(j1-j0-1)+len0)*sizeof(unsigned int));//allocate space for number of new edges in the resulting (after node delteion) mesh and all the remaining nodes;
          (*cl)=(unsigned int*)malloc(((j1-j0)*(j1-j0-1)+len0)*sizeof(unsigned int));
          (*vl)=(double*)malloc(((j1-j0)*(j1-j0-1)+len0)*sizeof(double));
          for(;k0<len0 && k1<cnt_curr;){
          	if(k0==j0) k0=j1;
          	if(rw0[k0]==nds_td0[j]){k0++; continue;}
          	if((cl0[k0]<cl_msh[k1]) || (cl0[k0]==cl_msh[k1] && rw0[k0]<rw_msh[k1])){
          		(*cl)[cnt_m]=cl0[k0];
          		(*rw)[cnt_m]=rw0[k0];
          		(*vl)[cnt_m]=vl0[k0];
          		cnt_m++; k0++;
          		
          	
          	}
          	else if(cl0[k0]==cl_msh[k1] && rw0[k0]==rw_msh[k1]){
          		(*cl)[cnt_m]=cl0[k0];
          		(*rw)[cnt_m]=rw0[k0];
          		(*vl)[cnt_m]=vl0[k0]+vl_msh[k1];
          		cnt_m++; k0++; k1++; 
          	
          	}
          	else{
          		(*cl)[cnt_m]=cl_msh[k1];
          		(*rw)[cnt_m]=rw_msh[k1];
          		(*vl)[cnt_m]=vl_msh[k1];
          		cnt_m++; k1++;
          	
          	}
          	
          
          }
          for(;k0<len0;){
          	if(k0==j0) k0=j1;
          	if(rw0[k0]==nds_td0[j]){k0++; continue;}
     		(*cl)[cnt_m]=cl0[k0];
     		(*rw)[cnt_m]=rw0[k0];
     		(*vl)[cnt_m]=vl0[k0];
     		cnt_m++; k0++;
          
          
          }
          for(;k1<cnt_curr;){
			(*cl)[cnt_m]=cl_msh[k1];
			(*rw)[cnt_m]=rw_msh[k1];
			(*vl)[cnt_m]=vl_msh[k1];
			cnt_m++; k1++;
          		
          
          }
          
          j++; i=0;
          //nds_i_val=i;
          
          tmp_p_ui=cl0; cl0=(*cl); free(tmp_p_ui);
          tmp_p_ui=rw0; rw0=(*rw); free(tmp_p_ui);
          tmp_p_d=vl0; vl0=(*vl); free(tmp_p_d);
          rw0=(unsigned int*)realloc(rw0,cnt_m*sizeof(unsigned int));
          cl0=(unsigned int*)realloc(cl0,cnt_m*sizeof(unsigned int));
          vl0=(double*)realloc(vl0,cnt_m*sizeof(double));
          len0=cnt_m;
          clock_gettime(CLOCK_MONOTONIC,&curr_time); tick1 = curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec;
          //fprintf(fp,"\ndt:,%llu",tick1-tick);
          *mode2+=tick1-tick;
          /*
          
          
          unsigned int k0=0,k1=0,k2=cnt_curr,k3=j1, cnt_m=0;//k0:[0;j0], k1:[0,cnt_curr], k2: [cnt_curr;2*cnt_curr], k3: [j1,len0];
          if(k0==j0) k0=j1;
          (*rw)=(unsigned int*)realloc((*rw),((len0)*((len0)-1))*2*sizeof(unsigned int));//allocate space for maximum number of new edges (if all nodes were connected);
          (*cl)=(unsigned int*)realloc((*cl),((len0)*((len0)-1))*2*sizeof(unsigned int));
          (*vl)=(double*)realloc((*vl),((len0)*((len0)-1))*2*sizeof(double));
          while(k1<cnt_curr && k2<2*cnt_curr && k0<len0){
               if(k0<len0 && rw0[k0]==nds_td0[j]) {k0++; continue;}
               if(cl_msh[k1]<cl_msh[k2]){
                    if(cl_msh[k1]<cl0[k0] || (cl_msh[k1]==cl0[k0] && rw_msh[k1]<rw0[k0])){
                         (*cl)[cnt_m]=cl_msh[k1];
                         (*rw)[cnt_m]=rw_msh[k1];
                         (*vl)[cnt_m]=vl_msh[k1];
                         cnt_m++;
                         k1++;
                         
                    }
                    else if(cl_msh[k1]==cl0[k0] && rw_msh[k1]==rw0[k0]){
                         (*cl)[cnt_m]=cl_msh[k1];
                         (*rw)[cnt_m]=rw_msh[k1];
                         (*vl)[cnt_m]=vl_msh[k1]+vl0[k0];
                         cnt_m++;
                         k1++;
                         k0++;
                         if(k0==j0) k0=j1;
                         
                         
                    }
                    else if((cl_msh[k1]==cl0[k0] && rw_msh[k1]>rw0[k0]) || cl_msh[k1]>cl0[k0]){
                         (*cl)[cnt_m]=cl0[k0];
                         (*rw)[cnt_m]=rw0[k0];
                         (*vl)[cnt_m]=vl0[k0];
                         cnt_m++;
                         k0++;
                         if(k0==j0) k0=j1;
                    
                    }
               
               }
               else if(cl_msh[k1]>=cl_msh[k2]){//part of array with k1 index has rows bigger than columns, part with k2 has columns bigger than rows;
                    if(cl_msh[k2]<cl0[k0] || (cl_msh[k2]==cl0[k0] && rw_msh[k2]<rw0[k0])){
                         (*cl)[cnt_m]=cl_msh[k2];
                         (*rw)[cnt_m]=rw_msh[k2];
                         (*vl)[cnt_m]=vl_msh[k2];
                         cnt_m++;
                         k2++;
                         
                    
                    }
                    else if(cl_msh[k2]==cl0[k0] && rw_msh[k2]==rw0[k0]){
                         (*cl)[cnt_m]=cl_msh[k2];
                         (*rw)[cnt_m]=rw_msh[k2];
                         (*vl)[cnt_m]=vl_msh[k2]+vl0[k0];
                         cnt_m++;
                         k2++;
                         k0++;
                         if(k0==j0) k0=j1;
                         
                         
                    
                    }
                    else if((cl_msh[k2]==cl0[k0] && rw_msh[k2]>rw0[k0]) || cl_msh[k2]>cl0[k0]){
                         (*cl)[cnt_m]=cl0[k0];
                         (*rw)[cnt_m]=rw0[k0];
                         (*vl)[cnt_m]=vl0[k0];
                         cnt_m++;
                         k0++;
                         if(k0==j0) k0=j1;
                    
                    }
               
               }
          
          }
          //////////////////////////////////////
          while(k1==cnt_curr && k2<2*cnt_curr && k0<len0){
               if(k0<len0 && rw0[k0]==nds_td0[j]) {k0++; continue;}
               if(cl_msh[k2]<cl0[k0] || (cl_msh[k2]==cl0[k0] && rw_msh[k2]<rw0[k0])){
                    (*cl)[cnt_m]=cl_msh[k2];
                    (*rw)[cnt_m]=rw_msh[k2];
                    (*vl)[cnt_m]=vl_msh[k2];
                    cnt_m++;
                    k2++;
                    
               }
               else if(cl_msh[k2]==cl0[k0] && rw_msh[k2]==rw0[k0]){
                    (*cl)[cnt_m]=cl_msh[k2];
                    (*rw)[cnt_m]=rw_msh[k2];
                    (*vl)[cnt_m]=vl_msh[k2]+vl0[k0];
                    cnt_m++;
                    k2++;
                    k0++;
                    if(k0==j0) k0=j1;
                    
               
               }
               else if((cl_msh[k2]==cl0[k0] && rw_msh[k2]>rw0[k0]) || cl_msh[k2]>cl0[k0]){
                    (*cl)[cnt_m]=cl0[k0];
                    (*rw)[cnt_m]=rw0[k0];
                    (*vl)[cnt_m]=vl0[k0];
                    cnt_m++;
                    k0++;
                    if(k0==j0) k0=j1;
               
               }
          }
          
          while(k1<cnt_curr && k2==2*cnt_curr && k0<len0){
               if(k0<len0 && rw0[k0]==nds_td0[j]) {k0++; continue;}
               if(cl_msh[k1]<cl0[k0] || (cl_msh[k1]==cl0[k0] && rw_msh[k1]<rw0[k0])){
                    (*cl)[cnt_m]=cl_msh[k1];
                    (*rw)[cnt_m]=rw_msh[k1];
                    (*vl)[cnt_m]=vl_msh[k1];
                    cnt_m++;
                    k1++;
                    
               
               }
               else if(cl_msh[k1]==cl0[k0] && rw_msh[k1]==rw0[k0]){
                    (*cl)[cnt_m]=cl_msh[k1];
                    (*rw)[cnt_m]=rw_msh[k1];
                    (*vl)[cnt_m]=vl_msh[k1]+vl0[k0];
                    cnt_m++;
                    k1++;
                    k0++;
                    if(k0==j0) k0=j1;
                    
               
               }
               else if((cl_msh[k1]==cl0[k0] && rw_msh[k1]>rw0[k0]) || cl_msh[k1]>cl0[k0]){
                    (*cl)[cnt_m]=cl0[k0];
                    (*rw)[cnt_m]=rw0[k0];
                    (*vl)[cnt_m]=vl0[k0];
                    cnt_m++;
                    k0++;
                    if(k0==j0) k0=j1;
               
               }
          
          }
          
          while(k1<cnt_curr && k2<2*cnt_curr && k0==len0){
               if(cl_msh[k1]<cl_msh[k2]){
                    (*cl)[cnt_m]=cl_msh[k1];
                    (*rw)[cnt_m]=rw_msh[k1];
                    (*vl)[cnt_m]=vl_msh[k1];
                    cnt_m++;
                    k1++;
               
               }
               else{//cl_msh[k1]==cl_msh[k2] or cl_msh[k1]>cl_msh[k2]; part of array with k1 index has rows bigger than columns, part with k2 has columns bigger than rows;
                    (*cl)[cnt_m]=cl_msh[k2];
                    (*rw)[cnt_m]=rw_msh[k2];
                    (*vl)[cnt_m]=vl_msh[k2];
                    cnt_m++;
                    k2++;
                    
               
               }
               
          
          }
          //////////////////////////
          //if((cnt_m+len0)>len0){
               rw0=(unsigned int*)realloc(rw0,2*(cnt_m+len0)*sizeof(unsigned int));
               cl0=(unsigned int*)realloc(cl0,2*(cnt_m+len0)*sizeof(unsigned int));
               vl0=(double*)realloc(vl0,2*(cnt_m+len0)*sizeof(double));
          
          //}
          while(k1==cnt_curr && k2==2*cnt_curr && k0<len0){
               if(rw0[k0]==nds_td0[j]) {k0++; continue;}
               (*cl)[cnt_m]=cl0[k0];
               (*rw)[cnt_m]=rw0[k0];
               (*vl)[cnt_m]=vl0[k0];
               cnt_m++;
               k0++;
               if(k0==j0) k0=j1;
     
          }
          while(k1<cnt_curr && k2==2*cnt_curr && k0==len0){
               (*cl)[cnt_m]=cl_msh[k1];
               (*rw)[cnt_m]=rw_msh[k1];
               (*vl)[cnt_m]=vl_msh[k1];
               cnt_m++;
               k1++;
          
          }
          while(k1==cnt_curr && k2<2*cnt_curr && k0==len0){
               (*cl)[cnt_m]=cl_msh[k2];
               (*rw)[cnt_m]=rw_msh[k2];
               (*vl)[cnt_m]=vl_msh[k2];
               cnt_m++;
               k2++;
          
          }
          ///////////////////
          j++;
          //nds_i_val=i;
          
          tmp_p_ui=(*cl); (*cl)=cl0; cl0=tmp_p_ui;// free(*cl);
          tmp_p_ui=(*rw); (*rw)=rw0; rw0=tmp_p_ui;// free(*rw);
          tmp_p_d=(*vl); (*vl)=vl0; vl0=tmp_p_d;// free(*vl);
          
          if(cnt_m>len0){
               rw0=(unsigned int*)realloc(rw0,2*cnt_m*sizeof(unsigned int));
               cl0=(unsigned int*)realloc(cl0,2*cnt_m*sizeof(unsigned int));
               vl0=(double*)realloc(vl0,2*cnt_m*sizeof(double));
          
          }
          len0=cnt_m;
          i=0;
          clock_gettime(CLOCK_MONOTONIC,&curr_time); tick1 = curr_time.tv_sec * 1000000000ll + curr_time.tv_nsec;
          //fprintf(fp,"\ndt:,%llu",tick1-tick);
          *mode2+=tick1-tick;*/
     
     
     }//////////////////////////////////////
     
     //free(*cl); free(*rw); free(*vl);
     (*cl)=cl0; (*rw)=rw0; (*vl)=vl0;
     //(*rw)=(unsigned int*)realloc((*rw),cnt_m*sizeof(unsigned int));
     //(*cl)=(unsigned int*)realloc((*cl),cnt_m*sizeof(unsigned int));
     //(*vl)=(double*)realloc((*vl),cnt_m*sizeof(double));
     *ln=cnt_m;
    
     
     free(nds_td0); //free(nds_td1); //free(rw0); free(cl0); free(vl0);
     free(rw_msh); free(cl_msh); free(vl_msh);
 
     //fclose(fp);


}

/*void group_deduct(unsigned int *row, unsigned int ***rw, unsigned int *col, unsigned int ***cl, double *val, double ***vl, unsigned int *len, unsigned int **ln, unsigned int *nds_td, unsigned int *nds_n, unsigned int* 

void sys_reduct(unsigned int *row,unsigned int** rw, unsigned int *col,unsigned int** cl, double *val, double** vl, unsigned int *len,unsigned int *ln, unsigned int *nds_td, unsigned int *nds_n, unsigned int* sz){TODO mode3;
     unsigned int **nds_grp;
     nds_grp=(unsigned int**) malloc((nds_n/(*sz)*sizeof(unsigned int*));
     unsigned int cnt_tmp=0;
     for(unsigned int i=0;i<(nds_n/(*sz);i++){
          nds_grp[i]=(unsigned int*)malloc((*sz)*sizeof(unsigned int));
          for(unsigned int j=0;j<*sz;j++){
               nds_grp[i][j]=nds_td[cnt_tmp];
               cnt_tmp++;
          
          }
     
     }
     #pragma omp parallel for num_threads((*nds_n)<100?(*nds_n):100)
     mode_2_del(row, rw, col, cl,val, vl, len,ln, unds_grp[i], *sz);
     
     

}*/

