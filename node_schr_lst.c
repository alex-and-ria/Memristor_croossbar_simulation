
































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
     
     unsigned int** npt=(unsigned int **) malloc(cl0[(*len)-1]*sizeof(unsigned int*));//table (array) that stores pointers to structs that correspond to nodes, node number i should have addres of the struct stored at npt[i];
     for(unsigned int i=0;i<cl0[*len-1];i++){
          npt[i]=NULL;
     
     }
     for(unsigned int i=0;i<*len;i++){//TODO form struct, record npt;
     
     }
     
     
     
     free(nds_td0); free(nds_td1); free(rw0); free(cl0); free(vl0);
     






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

