
































#include<stdio.h>


typedef struct node{
     unsigned int num;
     unsigned int n;
     unsigned int* nums;
     struct node* next;
     
} node;

#include"mode_2_3_alg.c"

void mode_1_alg(unsigned int *row,unsigned int** rw, unsigned int *col,unsigned int** cl, double *val,double** vl, unsigned int len,unsigned int *ln, unsigned int *nds_td, unsigned int nds_n, double thr_koef, unsigned char out_fl){
     unsigned int max_nds=col[len-1];
     //struct edge;
     /*typedef struct node{
          unsigned int num;
          unsigned int n;
          unsigned int* nums;
          struct node* next;
          
     } node;*/
     /*typedef struct edge{
          node* n1;
          node* n2;
          double val;
     
     } edge;*/
     node** node_arr=(node**) malloc(max_nds*sizeof(node*));
     //edge** edge_arr=(edge**) malloc(max_nds*max_nds*sizeof(edge*));
     double* edge_arr=(double*) malloc(max_nds*max_nds*sizeof(double));
     for(unsigned int i=0;i<max_nds*max_nds;i++){
          edge_arr[i]=-1;
     
     }
     unsigned int* nds_td0=(unsigned int*) malloc(max_nds*sizeof(unsigned int));
     unsigned int* nds_td_rem=(unsigned int*) malloc(max_nds*sizeof(unsigned int));
     unsigned int* nds_td2=(unsigned int*) malloc(max_nds*sizeof(unsigned int));
     unsigned int *buff_num=(unsigned int*) malloc(max_nds*sizeof(unsigned int));
     unsigned int *buff_mrg=(unsigned int*) malloc(max_nds*sizeof(unsigned int));
     double *sums=(double*)malloc(nds_n*sizeof(double));
      
     node* node_hd=(node*) malloc(1*sizeof(node));
     node_hd->num=col[0]; node_hd->nums=(unsigned int*)malloc((max_nds-1)*sizeof(unsigned int));
     node* curr_node=node_hd;
     node_arr[0]=curr_node;
     for(unsigned int i=col[0]+1;i<=col[len-1];i++){//assumption here is that node numbering is sequential without skipping the numbers;
          node_arr[i-1]=(node*)malloc(1*sizeof(node));
          node_arr[i-1]->num=i;
          node_arr[i-1]->nums=(unsigned int*)malloc((max_nds-1)*sizeof(unsigned int));
          node_arr[i-2]->next=node_arr[i-1];
          
     }
     //curr_node->next=NULL;
     node_arr[col[len-1]-1]->next=NULL;
     ///////////////////////////////at this point nodes should be set up;
     
     
     
     
     //edge *curr_edge=NULL;
     unsigned int curr_col=col[0];
     unsigned int i0=0;
     for(unsigned int i=1;i<len;i++){
          if(curr_col!=col[i]){
               node_arr[col[i-1]-1]->n=i-i0;
               for(unsigned int j=i0;j<i;j++){
                    node_arr[col[i-1]-1]->nums[j-i0]=row[j];
                    /*curr_edge=(edge*)malloc(1*sizeof(edge));
                    if(col[i-1]<row[j]){
                         curr_edge->n1=node_arr[col[i-1]-1];
                         curr_edge->n2=node_arr[row[j]-1];
                         
                    }
                    else{
                         curr_edge->n1=node_arr[row[j]-1];
                         curr_edge->n2=node_arr[col[i-1]-1];
                    
                    }
                    curr_edge->val=val[j];*/
                    edge_arr[(col[i-1]-1)*max_nds+row[j]-1]=val[j];//curr_edge;
                    edge_arr[(row[j]-1)*max_nds+(col[i-1]-1)]=val[j];//curr_edge;
               
               }
               i0=i;
               curr_col=col[i];
          }
          
     
     }
     node_arr[max_nds-1]->n=len-i0;
     for(unsigned int j=i0;j<len;j++){
          node_arr[max_nds-1]->nums[j-i0]=row[j];
          /*curr_edge=(edge*)malloc(1*sizeof(edge));
          curr_edge->n1=node_arr[row[j]-1];
          curr_edge->n2=node_arr[max_nds-1];
          curr_edge->val=val[j];*/
          edge_arr[(max_nds-1)*max_nds+row[j]-1]=val[j];//curr_edge;
          edge_arr[(row[j]-1)*max_nds+(max_nds-1)]=val[j];//curr_edge;
     
     }
     /////////////////////////////////////at this point edges should be set up too;
     
     
     
     
     //unsigned int* nds_td0=(unsigned int*) malloc(max_nds*sizeof(unsigned int));
     
     unsigned int nds_n0=0,nds_n_rem=0,nds_n2=nds_n;
     for(unsigned int i=0;i<nds_n;i++){
          nds_td2[i]=nds_td[i];
     
     }
     
     //unsigned int dbg_cnt=0,dbg_max=4;
     double curr_thr_koef=0;
     while(1/*dbg_cnt<dbg_max*/){
     for(unsigned int i=0;i<max_nds;i++){
          nds_td0[i]=0; nds_td_rem[i]=0;
     }
     for(unsigned int i=0;i<nds_n2;i++){
          nds_td0[nds_td2[i]-1]=nds_td2[i]; nds_td_rem[nds_td2[i]-1]=nds_td2[i];
     }
     /*for(curr_node=node_hd;curr_node!=NULL;curr_node=curr_node->next){
          if(nds_td0[curr_node->num-1]!=0){
               for(unsigned int i=0;i<curr_node->n;i++){
                    if(curr_node->num < curr_node->nums[i]){
                         nds_td0[curr_node->nums[i]-1]=0;
                    
                    }
                    for(j=0;j<node_arr[curr_node->nums[i]-1]->n;j++){
                         if(curr_node->num < node_arr[curr_node->nums[i]-1]->nums[j]){
                              nds_td0[node_arr[curr_node->nums[i]-1]->nums[j]-1]=0;
                         
                         }
                         
                    
                    }
               
               }
          
          }
          
     
     }*/
     for(unsigned int i=0;i<max_nds;){
          while(i<max_nds && nds_td0[i]==0){
               i++;
          
          }
          if(i<max_nds){
               curr_node=node_arr[i];//nds_td0[i] should be equal to i+1 if it is not zero; curr_node->num should also be equal to i+1;
               for(unsigned int j=0;j<curr_node->n;j++){
                    if(nds_td0[curr_node->nums[j]-1]!=0 && nds_td0[i] < curr_node->nums[j]) nds_td0[curr_node->nums[j]-1]=0;
                    for(unsigned int k=0;k<node_arr[curr_node->nums[j]-1]->n;k++){
                         if(nds_td0[node_arr[curr_node->nums[j]-1]->nums[k]-1]!=0 && nds_td0[i] < node_arr[curr_node->nums[j]-1]->nums[k]){
                              nds_td0[node_arr[curr_node->nums[j]-1]->nums[k]-1]=0;
                         
                         }
                    
                    }
               
               }
               i++;
          }
          
     }
     for(unsigned int i=0;i<max_nds;i++){
          if(nds_td0[i]!=0){
               nds_td_rem[nds_td0[i]-1]=0;
               nds_td0[nds_n0]=nds_td0[i];
               nds_n0++;
          
          }
          if(nds_td_rem[i]!=0){
               nds_td_rem[nds_n_rem]=nds_td_rem[i];
               nds_n_rem++;
          
          }
     
     }
     if(nds_n0>=1) curr_thr_koef=(nds_n0+0.)/(nds_n0+nds_n_rem);
     if(nds_n0<=1 || curr_thr_koef<thr_koef){
          break;
     
     }
     /////////////////////////////////////////at this point nds_td0 (nodes to delete) and nds_td_rem (nodes left, delete on next iteration) should set up;
     
     
     
     
     unsigned int edge_idx0=0,edge_idx1=0,edge_idx2=0,n1,n2;
     unsigned int bf_nm_cnt=0,bf_mg_cnt=0;
     unsigned int *ui_ptr;
    
     for(unsigned int i=0;i<nds_n0;i++){
          sums[i]=0;
          n1=nds_td0[i];
          for(unsigned int j=0;j<node_arr[n1-1]->n;j++){
               n2=node_arr[n1-1]->nums[j];
               sums[i]+=edge_arr[(n1-1)*max_nds+(n2-1)];//->val;
               
          
          }
     }
     
     for(unsigned int i=0;i<nds_n0;i++){
          for(unsigned int j=0;j<node_arr[nds_td0[i]-1]->n;j++){
               n1=node_arr[nds_td0[i]-1]->nums[j];
               bf_nm_cnt=0;
               for(unsigned int k=0;k<node_arr[nds_td0[i]-1]->n;k++){
                    if(j==k) continue;
                    n2=node_arr[nds_td0[i]-1]->nums[k];
                    buff_num[bf_nm_cnt]=n2; bf_nm_cnt++;
                    edge_idx0=(n1-1)*max_nds+n2-1;
                    edge_idx1=(node_arr[nds_td0[i]-1]->num-1)*max_nds+n1-1;
                    edge_idx2=(node_arr[nds_td0[i]-1]->num-1)*max_nds+n2-1;
                    //edge_idx1=node_arr[nds_td0[i]-1]->nums[k]-1)*max_nds+node_arr[nds_td0[i]-1]->nums[j]-1;
                    /*curr_idx0=(node_arr[nds_td0[i]-1]->n-1)*j+k-1;
                    curr_idx1=(node_arr[nds_td0[i]-1]->n-1)*k+j;
                    rw_buff[curr_idx0]=node_arr[nds_td0[i]-1]->nums[j];
                    cl_buff[curr_idx0]=node_arr[nds_td0[i]-1]->nums[k];
                    rw_buff[curr_idx1]=node_arr[nds_td0[i]-1]->nums[k];
                    cl_buff[curr_idx1]=node_arr[nds_td0[i]-1]->nums[j];*/
                    if(n1<n2 && edge_arr[edge_idx0]>0/*!=NULL*/){//this edge between node nums[j] and nums[k];//if there are no common neighbour between two nodes, their neighbours do not form common edge, but there might be this edge in original graph;
                         //edge_arr[edge_idx0]->val+=(edge_arr[edge_idx1]->val*edge_arr[edge_idx2]->val)/(sums[i]);
                         edge_arr[edge_idx0]+=(edge_arr[edge_idx1]*edge_arr[edge_idx2])/(sums[i]);
                         //edge_arr[edge_idx1]=edge_arr[edge_idx0];
                    }
                    else if(n1<n2){
                         /*curr_edge=(edge*)malloc(1*sizeof(edge));
                         curr_edge->n1=node_arr[n1-1]; curr_edge->n2=node_arr[n2-1];
                         curr_edge->val=(edge_arr[edge_idx1]->val*edge_arr[edge_idx2]->val)/(sums[i]);*/
                         edge_arr[edge_idx0]=(edge_arr[edge_idx1]*edge_arr[edge_idx2])/(sums[i]);
                         //edge_arr[edge_idx0]=curr_edge;//case when n1 and n2 swap, would be one of the iterations too;
                    
                    }
                    else{
                         edge_arr[edge_idx0]=edge_arr[(n2-1)*max_nds+n1-1];
                    
                    }
               }
               bf_mg_cnt=0;
               unsigned int ll=0,qq=0;
               for(; ll<node_arr[n1-1]->n && qq<bf_nm_cnt;){
                    if(node_arr[n1-1]->nums[ll] <= buff_num[qq]){
                         if(node_arr[n1-1]->nums[ll]==nds_td0[i]){//for node node_td0[i] and its neighbours nums=node_arr[node_td[i]-1]->nums, each of nodes in nums has node_td0[i] as their own neighbour (if a has edge with b then b has edge with a), so skip nds_td0[i] node in its neighbour list of neighbours; no other node from nds_td0 is in neigbour list because othewise this node whould be the common node (common neigbour) for two nodes from nds_td0;
                              ll++;
                         }
                         else{
                              buff_mrg[bf_mg_cnt]=node_arr[n1-1]->nums[ll]; bf_mg_cnt++;
                              if(node_arr[n1-1]->nums[ll] == buff_num[qq]) qq++;
                              ll++;
                         }
                    }
                    else{//if node_arr[n1-1]->nums[ll] == buff_num[qq] it means that this is old node from n1 to n2 that is in parralel to new one; parrallel conductances were summed; update: parralel conductance does not appear if no common nodes;
                         buff_mrg[bf_mg_cnt]=buff_num[qq]; bf_mg_cnt++; qq++;
                         
                    }
               
               }
               while(ll<node_arr[n1-1]->n){
                    if(node_arr[n1-1]->nums[ll]==nds_td0[i]){
                         ll++;
                    }
                    else{
                         buff_mrg[bf_mg_cnt]=node_arr[n1-1]->nums[ll]; bf_mg_cnt++; ll++;
                    }
               
               }
               while(qq<bf_nm_cnt){
                    buff_mrg[bf_mg_cnt]=buff_num[qq]; bf_mg_cnt++; qq++;
               
               }
               ui_ptr=node_arr[n1-1]->nums; node_arr[n1-1]->nums=buff_mrg; buff_mrg=ui_ptr;//intergrity is supported by nums arrays;
               node_arr[n1-1]->n=bf_mg_cnt;
          
          }
          /*n1=nds_td0[i];
          for(unsigned int j=0;j<node_arr[nds_td0[i]-1]->n;j++){
               n2=node_arr[nds_td0[i]-1]->nums[j];
               //free(edge_arr[(n1-1)*max_nds+n2-1]);
               edge_arr[(n1-1)*max_nds+n2-1]=-1;//NULL;
               edge_arr[(n2-1)*max_nds+n1-1]=-1;//NULL;
          
          }*/
     
     }
     /////////////////////////////////////////////////////////////////////at this point star mesh transform should be completed for this iteration
     
     
     
     node* prev_node=NULL;curr_node=node_hd;
     for(unsigned int i=0; i<nds_n0 && curr_node!=NULL;){
          if(curr_node->num==nds_td0[i]){
               free(curr_node->nums);
               node_arr[curr_node->num-1]=NULL;
               if(prev_node==NULL){
                    node_hd=curr_node->next;
                    free(curr_node);
                    curr_node=node_hd;
               
               }
               else{
                    prev_node->next=curr_node->next;
                    free(curr_node);
                    curr_node=prev_node->next;
               
               }
               i++;
          
          }
          else if(curr_node->num<nds_td0[i]){
               /*unsigned int idxl=0, ordr=0;
               unsigned int j=0,k=0;
               for(j=0,k=0;j<curr_node->n && k<nds_n0;){
                    if(curr_node->nums[j] == nds_td0[k] && idxl==0){
                         idxl=j+1; j++; k++; ordr=1;
                    
                    }
                    else if(curr_node->nums[j] == nds_td0[k] && idxl>0){
                         for(unsigned int ll=idxl;ll<j;ll++){
                              curr_node->nums[ll-ordr]=curr_node->nums[ll];
                              idxl=j+1; ordr++;
                         
                         }
                         j++; k++;
                    
                    }
                    else if(curr_node->nums[j] < nds_td0[k]){
                         j++;
                    
                    }
                    else{
                         k++;
                    
                    }
                    
               
               }
               while(j<curr_node->n){
                    if(idxl==0){
                         break;//array is already fine, no nodes from nds_td0 found;
                    
                    }
                    else{
                         curr_node->nums[j-ordr]=curr_node->nums[j];
                         j++;
                         
                    }
               
               }
               //even if k<nds_n0, nodes from the nums array ended, hence they are fine;
               curr_node->n-=ordr;*/
               prev_node=curr_node;
               curr_node=curr_node->next;
          
          }
     }
     /*while(curr_node!=NULL){
          unsigned int idxl=0,ordr=0; unsigned int j=0,k=0;
          for(j=0,k=0;j<curr_node->n && k<nds_n0;){
               if(curr_node->nums[j] == nds_td0[k] && idxl==0){
                    idxl=j+1; j++; k++; ordr=1;
               }
               else if(curr_node->nums[j] == nds_td0[k] && idxl>0){
                    for(unsigned int ll=idxl;ll<j;ll++){
                         curr_node->nums[ll-ordr]=curr_node->nums[ll];
                         idxl=j+1; ordr++;
                    }
                    j++; k++;
               }
               else if(curr_node->nums[j] < nds_td0[k]){
                    j++;
               }
               else{
                    k++;
               }
          }
          while(j<curr_node->n){
               if(idxl==0){
                    break;
               }
               else{
                    curr_node->nums[j-ordr]=curr_node->nums[j];
                    j++;
               }
          }
          curr_node->n-=ordr;
          curr_node=curr_node->next;
     
     }*/
     /////////////////////////////////at this point nums array are edited to match new graph connectivity;
     
     
     
     ui_ptr=nds_td2; nds_td2=nds_td_rem; nds_td_rem=ui_ptr;
     nds_n2=nds_n_rem; nds_n0=0; nds_n_rem=0;
     //dbg_cnt++;
     }
     
     if(out_fl==1){
          *ln=0;
          for(curr_node=node_hd;curr_node!=NULL;curr_node=curr_node->next){
               *ln+=curr_node->n;
          
          }
          (*cl)=(unsigned int*) malloc((*ln)*sizeof(unsigned int));
          (*rw)=(unsigned int*) malloc((*ln)*sizeof(unsigned int));
          (*vl)=(double*) malloc((*ln)*sizeof(double));
          unsigned int cnt_out=0;
          for(curr_node=node_hd;curr_node!=NULL;curr_node=curr_node->next){
               for(unsigned int i=0;i<curr_node->n;i++){
                    (*cl)[cnt_out]=curr_node->num;
                    (*rw)[cnt_out]=curr_node->nums[i];
                    (*vl)[cnt_out]=edge_arr[(curr_node->num-1)*max_nds+curr_node->nums[i]-1];//->val;
                    cnt_out++;
                    
               
               }
          
          }
          free(nds_td0);
          free(nds_td_rem);
          free(sums);
          node* prev_node=NULL;
          for(curr_node=node_hd;curr_node!=NULL;){
               prev_node=curr_node;
               curr_node=curr_node->next;
               free(prev_node->nums);
               free(prev_node);
               
            
          
          }
          free(node_arr);
          free(edge_arr);
          
          
     }
     else{
          out_fl--;
          mode_2_alg(node_arr,edge_arr,nds_td2,nds_n2, &node_hd,max_nds,rw,cl,vl,ln,out_fl);
          //mode_2_alg(node** node_arr,double* edge_arr,unsigned int *nds_td, unsigned int nds_n, node** node_hd_, unsigned int max_nds, unsigned int** rw, unsigned int** cl, double **vl, unsigned int *ln, unsigned char out_fl);
          free(nds_td0);
          free(nds_td_rem);
          free(sums);
          free(node_arr);
          free(edge_arr);
          
          
     
     
     }
     
     free(nds_td2);
     free(buff_num);
     free(buff_mrg);
     
     

}

