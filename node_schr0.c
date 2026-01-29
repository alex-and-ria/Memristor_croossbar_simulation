
































#define th_nb_koef (0.8)


void node_analyzer(unsigned int* row, unsigned int *col, unsigned int *len, unsigned int* nds_td, unsigned int* nds_n, double* nb_koef){//adj_m starts count from node 1 (Matlab indexing); it is represented in sparse Matlab format via [row,col,val]; col are sorted; nds_td should count from one (Matlab idexing), and be sorted; it is a list of nodes to delete; check nodes and delete one of them from nds_td if two are neighbours;

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
}

void dense_rdct(unsigned int *row,unsigned int** rw, unsigned int *col,unsigned int** cl, double *val, double** vl, unsigned int *len,unsigned int *ln, unsigned int *nds_td, unsigned int *nds_n){
     double nb_koef=0;
     unsigned int *nds_td0=(unsigned int*)malloc((*nds_n)*sizeof(unsigned int)); unsigned int nds_td_cnt=0;
     unsigned int *nds_td1=(unsigned int*)malloc((*nds_n)*sizeof(unsigned int));
     for(unsigned int i=0;i<*nds_n;i++){
          nds_td0[i]=nds_td[i];
     }
     node_analyzer(row,col,len,nds_td,nds_n, &nb_koef);
     unsigned int nds_n_nz=0;
     unsigned int swp_tmp=0; double swp_tmp_d=0;
     /*(*rw)=(unsigned int*)malloc((*len)*sizeof(unsigned int));
     (*cl)=(unsigned int*)malloc((*len)*sizeof(unsigned int));
     (*vl)=(double*)malloc((*len)*sizeof(double));*/
     unsigned int *rw0=(unsigned int*)malloc((*len)*sizeof(unsigned int));
     unsigned int *cl0=(unsigned int*)malloc((*len)*sizeof(unsigned int));
     double *vl0=(double*)malloc((*len)*sizeof(double));
     unsigned int *tmp_p_ui; double *tmp_p_d;
     for(unsigned int i=0;i<*len;i++){
          cl0[i]=col[i]; rw0[i]=row[i]; vl0[i]=val[i];
     
     }
     //rw0=row; cl0=col; vl0=val;
     while (nb_koef<th_nb_koef){
          nds_n_nz=0; nds_td_cnt=0;
          for(unsigned int i=0;i<*nds_n;){//delete zeros; keep nodes to delete on next iterations;
               while(i<*nds_n && nds_td[i]==0){
                    nds_td1[nds_td_cnt]=nds_td0[i];
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
          star_mesh_base(rw0,rw,cl0,cl,vl0,vl,len,ln,nds_td0,nds_n);
          
          for(unsigned int i=0;i<*ln;i++){//sort with selection sort (major is column sorting, secondary is row sorting); faster alghourithms might be equivalently used (like mege sort, or quick sort);
               for(unsigned int j=i+1;j<*ln;j++){
                    if((*cl)[j]<(*cl)[i] || ((*cl)[j]==(*cl)[i] && (*rw)[j]<((*rw)[i]))){
                         swp_tmp=(*cl)[i]; (*cl)[i]=(*cl)[j]; (*cl)[j]=swp_tmp;
                         swp_tmp=(*rw)[i]; (*rw)[i]=(*rw)[j]; (*rw)[j]=swp_tmp;
                         swp_tmp_d=(*vl)[i]; (*vl)[i]=(*vl)[j]; (*vl)[j]=swp_tmp_d;
                     
                    }
               
               }
          
          }//at this point *rw,*cl, and *vl should be sorted like standart Matlab column major format;
          
          tmp_p_ui=(*cl); (*cl)=cl0; cl0=tmp_p_ui; free(*cl);
          tmp_p_ui=(*rw); (*rw)=rw0; rw0=tmp_p_ui; free(*rw);
          tmp_p_d=(*vl); (*vl)=vl0; vl0=tmp_p_d; free(*vl);
          
          
          *len=*ln;
          *nds_n=nds_td_cnt;//since number of nodes to delete is decreesing at each new iteration, no need to realloc memory;
          for(unsigned int i=0;i<*nds_n;i++){
               nds_td0[i]=nds_td1[i]; nds_td[i]=nds_td1[i];
          }
          node_analyzer(rw0,cl0,len,nds_td,nds_n, &nb_koef);
          
     
     
     
     }
     
     //nds_td0 nds_n rw0 cl0 vl0 len
     (*rw)=(unsigned int*)malloc(((*len)*((*len)-1))/2*sizeof(unsigned int));//allocate space for maximum number of new edges (if all nodes were connected);
     (*cl)=(unsigned int*)malloc(((*len)*((*len)-1))/2*sizeof(unsigned int));
     (*vl)=(double*)malloc(((*len)*((*len)-1))/2*sizeof(double));
     unsigned int *nz_nd_a=(unsigned int*) malloc(2*((*len)-(*nds_n))*sizeof(unsigned int));//approximate number that is sufficient to keep reamaining nodes;
     double cnds_sum=0; unsigned int cnt_curr=0,cnt_curr_init=0; unsigned int cnt_nz_nds=0; unsigned char fl_nz=0; unsigned int nds_i_val=0;
     for(unsigned int i=0,j=0,j0=0,j1=0;i<*len && j<*nds_n;){//*nds_n is smaller or equal to *len;
          while(cl0[i]<nds_td0[j]){
               if(fl_nz==0){
                    nz_nd_a[2*cnt_nz_nds]=i;
                    fl_nz=1;
               
               }
               i++;
          
          }
          if(fl_nz==1){
               fl_nz=0; nz_nd_a[2*cnt_nz_nds+1]=i; cnt_nz_nds++;//TODO!;
          
          }//set indexes for nodes that should not be deleted at this point;
          j0=i;
          cnds_sum=0;
          while(i<*len && cl0[i]==nds_td0[j]){
               cl0[i]=0;
               cnds_sum+=vl0[i];
               i++;
               
          
          
          
          }
          j1=i;//index boundearies for current node to delete;
          //parallel
          for(unsigned int kk=j0;kk<j1-1;kk++){
               cnt_curr=cnt_curr_init;
               for(unsigned int kk0=j1-j0-1;kk0!=j1-kk-1;kk0--){
                    cnt_curr+=kk0;
               
               }//calculating index to avoid dependencies with previous iterations;
               
               for(unsigned int ll=kk+1;ll<j1;ll++){
                    (*cl)[cnt_curr+ll-kk-1]=rw0[kk];//(*cl) gets smaller index;
                    (*rw)[cnt_curr+ll-kk-1]=rw0[ll];//(*rw) gets bigger index;
                    (*vl)[cnt_curr+ll-kk-1]=1./((1./vl0[kk])*(1./vl0[ll])*cnds_sum);
                    
               
               }
               
          
          }
          cnt_curr_init=cnt_curr+1;
          j++;
          nds_i_val=i;
     
     
     }//////////////////////////////////////
     
     //TODO_ delete duplicates *cl, merge (duplicates + union); sort;
     for(unsigned int i=0;i<cnt_curr_init;i++){
          if((*cl)[i]==0) continue;
          for(unsigned int j=i+1;j<cnt_curr_init;j++){
               if((*cl)[i]==(*cl)[j] && (*rw)[i]==(*rw)[j]){
                    (*vl)[i]+=(*vl)[j]; (*cl)[j]=0;
                    
               
               }
               
          
          }//merging parrallel nodes within (*cl), (*rw)
          
          for(unsigned int j=0;j<cnt_nz_nds;j++){
               for(unsigned int kk=nz_nd_a[2*j];kk<nz_nd_a[2*j+1];kk++){
                    if(rw0[kk]<cl0[kk] && (*cl)[i]==rw0[kk] && (*rw)[i]==cl0[kk]){//rw0[j]<cl0[j] can evaluate to true if cl0[j]!=0 (indexes of rw0 are greater or equal to 1);
                         (*vl)[i]+=vl0[kk];
                         cl0[kk]=0;
                    
                    }
               
               
               
               }
               
          
          }//merging parallel nodes with remaining (original) nodes, whose node number is smaller that the biggest node number of *nds_td;
          
          for(unsigned int j=nds_i_val;j<*len;j++){
               if(rw0[j]<cl0[j] && (*cl)[i]==rw0[j] && (*rw)[i]==cl0[j]){//rw0[j]<cl0[j] can evaluate to true if cl0[j]!=0 (indexes of rw0 are greater or equal to 1);
                    (*vl)[i]+=vl0[j];
                    cl0[j]=0;
               
               }
          }
     
     }//at this point parallel conductances should be merged;
     ///////////////////////////////////////////////
     unsigned int cl_idx_cnt=0;
     for(unsigned int i=0; i<cnt_curr_init && (*cl)[i]!=0;i++,cl_idx_cnt++){
     }
     
     for(unsigned int i=cl_idx_cnt;i<cnt_curr_init;){
          while(i<cnt_curr_init && (*cl)[i]==0){
               i++;
          
          }
          while(i<cnt_curr_init && (*cl)[i]!=0){
               (*cl)[cl_idx_cnt]=(*cl)[i];
               (*rw)[cl_idx_cnt]=(*rw)[i];
               (*vl)[cl_idx_cnt]=(*vl)[i];
               i++; cl_idx_cnt++;
          
          }
     }//at this point zeros should be removed;
     //////////////////////////////////////
     for(unsigned int i=0;i<cl_idx_cnt;i++){
          (*cl)[cl_idx_cnt+i]=(*rw)[i]; (*rw)[cl_idx_cnt+i]=(*cl)[i]; (*vl)[cl_idx_cnt+i]=(*vl)[i];
     
     }//recording twice to comply with adjacency matrix format;
     cl_idx_cnt*=2;
     for(unsigned int i=0;i<cnt_nz_nds;i++){//record nodes that remain after deletion (original ones);
          for(unsigned int j=nz_nd_a[2*j];j<nz_nd_a[2*j+1];j++){
               if(rw0[j]<cl0[j]){//rw0[j]<cl0[j] can evaluate to true if cl0[j]!=0 (indexes of rw0 are greater or equal to 1) (this if condition avoids cl0[j]==0); this condition alows to avoid counting same edge twice;
                    (*cl)[cl_idx_cnt]=rw0[j]; (*rw)[cl_idx_cnt]=cl0[j]; (*vl)[cl_idx_cnt]=vl0[j];
                    cl_idx_cnt++;
                    (*cl)[cl_idx_cnt]=cl0[j]; (*rw)[cl_idx_cnt]=rw0[j]; (*vl)[cl_idx_cnt]=vl0[j];
                    cl_idx_cnt++;//recored twice to comply with adjacecy matrix format;
                    
               
               }
          
          
          
          }
          
     
     }
     
     for(unsigned int i=nds_i_val;i<*len;i++){//record nodes that are bigger then the biggest index *nds_td, and are not duplicates;
          if(rw0[i]<cl0[i] ){//rw0[j]<cl0[j] can evaluate to true if cl0[j]!=0 (indexes of rw0 are greater or equal to 1);
               (*cl)[cl_idx_cnt]=rw0[i]; (*rw)[cl_idx_cnt]=cl0[i]; (*vl)[cl_idx_cnt]=vl0[i];
               cl_idx_cnt++;
               (*cl)[cl_idx_cnt]=cl0[i]; (*rw)[cl_idx_cnt]=rw0[i]; (*vl)[cl_idx_cnt]=vl0[i];
               cl_idx_cnt++;//recored twice to comply with adjacecy matrix format;
          
          }
     }
     ///////////////////////////////////////
     
     (*rw)=(unsigned int*)realloc((*rw),cl_idx_cnt*sizeof(unsigned int));
     (*cl)=(unsigned int*)realloc((*cl),cl_idx_cnt*sizeof(unsigned int));
     (*vl)=(double*)realloc((*vl),cl_idx_cnt*sizeof(double));
     *len=cl_idx_cnt;
     for(unsigned int i=0;i<*len;i++){//sort with selection sort (major is column sorting, secondary is row sorting); faster alghourithms might be equivalently used (like mege sort, or quick sort);
          for(unsigned int j=i+1;j<*len;j++){
               if((*cl)[j]<(*cl)[i] || ((*cl)[j]==(*cl)[i] && (*rw)[j]<((*rw)[i]))){
                    swp_tmp=(*cl)[i]; (*cl)[i]=(*cl)[j]; (*cl)[j]=swp_tmp;
                    swp_tmp=(*rw)[i]; (*rw)[i]=(*rw)[j]; (*rw)[j]=swp_tmp;
                    swp_tmp_d=(*vl)[i]; (*vl)[i]=(*vl)[j]; (*vl)[j]=swp_tmp_d;
                
               }
          
          }
     
     }//at this point *rw,*cl, and *vl should be sorted like standart Matlab column major format;
     
     
     
     
     
     
     free(nds_td0); free(nds_td1); free(rw0); free(cl0); free(vl0); free(nz_nd_a);
     






}

/*void sys_reduct(unsigned int *row,unsigned int** rw, unsigned int *col,unsigned int** cl, double *val, double** vl, unsigned int *len,unsigned int *ln, unsigned int *nds_td, unsigned int *nds_n, unsigned int* sz){
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

