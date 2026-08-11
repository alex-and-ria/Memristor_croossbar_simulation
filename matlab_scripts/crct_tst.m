

































m=256; n=m; batch_size=2;
Gwl=1./100; Gbl=4./100;
[G_adj, Vin, Cnds]=init_cb(m,n,batch_size,Gwl,Gbl,0);
[row,col,val]=find(G_adj);
[G_m,Ivec]=adj_to_lapl(G_adj,m,n,Vin);
     [L,U,P]=lu(G_m);
     y=L\(P*Ivec);
     x_orgn=U\y;

%figure(1); spy(G_m);
%figure(2); spy(L);
%figure(3); spy(U);

loadlibrary('../libnode_schr.so','../node_schr.h')
%libfunctions('libnode_schr','-full')
nds_td=1:2*m*n; [nds_td,nds_tgt]=pune_ntd(nds_td,m,n);

row_p=libpointer('uint32Ptr',row);
rw_vp=libpointer('uint64Ptr',0);%memory for 64-bit address, to keep raw address (for tripple pointer);
col_p=libpointer('uint32Ptr',col);
cl_vp=libpointer('uint64Ptr',0);
val_p=libpointer('doublePtr',val);
vl_vp=libpointer('uint64Ptr',0);
len_p=libpointer('uint32Ptr',length(row));
len_pp=libpointer('uint32PtrPtr');
nds_td_p=libpointer('uint32Ptr',nds_td); nds_n=libpointer('uint32Ptr',length(nds_td)); 
nds_td1_p=libpointer('uint32Ptr',nds_tgt); nds_n1=length(nds_tgt); 
n_th_p=libpointer('uint32Ptr',0);

%null_uip=libpointer('uint32Ptr');
%max_m_sz=nds_n1;%set 2 submatrixes for mode 3;


nds_n_req=length(nds_td)+3;%nds_n_req<(length(nds_td)+length(nds_tgt))
if(nds_n_req>length(nds_td))
     max_m_sz=length(nds_tgt)-(nds_n_req-length(nds_td));
else
     max_m_sz=length(nds_tgt);
end
calllib('libnode_schr','dense_rdct',row_p,rw_vp,...
     col_p,cl_vp,...
     val_p,vl_vp,...
     len_p,len_pp,...
     nds_td_p,nds_n,...
     0.,...%smaller the number for mode_1_2_3, the more mode 1 is utilized;
     nds_td1_p, nds_n1,...
     n_th_p,max_m_sz,...%set 2 submatrixes for mode 3;
     2,0,...;%mode_1_2_3
     m,n,nds_n_req)

setdatatype(len_pp.Value,'uint32Ptr',n_th_p.Value);
if(nds_n_req<=length(nds_td))%means that only mode 1 and, maybe mode 2 is called;
     rw_c=libpointer('uint32PtrPtr');
     cl_c=libpointer('uint32PtrPtr');
     vl_c=libpointer('doublePtrPtr');
     calllib('libnode_schr','get_dbg_arr', rw_vp, cl_vp, vl_vp,...
               rw_c,cl_c,vl_c,0);%store data in rw_c,cl_c,vl_c;
     setdatatype(rw_c.Value,'uint32Ptr',len_pp.Value(1));
     setdatatype(cl_c.Value,'uint32Ptr',len_pp.Value(1));
     setdatatype(vl_c.Value,'doublePtr',len_pp.Value(1));
     G_iter=sparse(rw_c.Value,cl_c.Value,vl_c.Value);
     [G_m, Ivec0]=adj_to_lapl(G_iter,m,n,Vin);
     [L,U,P]=lu(G_m);
     y=L\(P*Ivec0);
     x_curr=U\y;
     sol_diff=x_orgn-x_curr;
     max(abs(sol_diff(~isnan(sol_diff))))

     
elseif(nds_n_req>length(nds_td))
     for(ii=0:(n_th_p.Value-1))
          rw_c=libpointer('uint32PtrPtr');
          cl_c=libpointer('uint32PtrPtr');
          vl_c=libpointer('doublePtrPtr');
          calllib('libnode_schr','get_dbg_arr', rw_vp, cl_vp, vl_vp,...
               rw_c,cl_c,vl_c,ii);%store data in rw_c,cl_c,vl_c;
          setdatatype(rw_c.Value,'uint32Ptr',len_pp.Value(ii+1));
          setdatatype(cl_c.Value,'uint32Ptr',len_pp.Value(ii+1));
          setdatatype(vl_c.Value,'doublePtr',len_pp.Value(ii+1));
          G_iter=sparse(rw_c.Value,cl_c.Value,vl_c.Value);
          [G_m, Ivec0]=adj_to_lapl(G_iter,m,n,Vin);
          [L,U,P]=lu(G_m); y=L\(P*Ivec0); x_curr=U\y;
          sol_diff=x_orgn-x_curr;
          max(abs(sol_diff(~isnan(sol_diff))))
          
     end
     
     
     
end


calllib('libnode_schr','data_free',rw_vp, cl_vp, vl_vp,len_pp,n_th_p);
%q=33;
unloadlibrary libnode_schr














function G_adj=adj_m_w(Cnds,Gwl,Gbl)
     %node numberring is as n_nd puts them; then nodes with soudetrces follow,
     %then ground;
     m=size(Cnds,1); n=size(Cnds,2);
     G_adj=sparse(1,1,0,2*m*n+m+1,2*m*n+m+1,m*(5+3*(n-2))+n*(5+3*(m-2))+m+1);
     n_nd(0,0,0,n);
     for ii=1:m
          for jj=1:n
               if(jj==1)
                    G_adj(n_nd(ii,jj,0),2*m*n+ii)=Gwl;
                    G_adj(n_nd(ii,jj,0),n_nd(ii,jj+1,0))=Gwl;
                    G_adj(n_nd(ii,jj,0),n_nd(ii,jj,1))=Cnds(ii,jj);
                    
               elseif(jj==n)
                    G_adj(n_nd(ii,jj,0),n_nd(ii,jj-1,0))=Gwl;
                    G_adj(n_nd(ii,jj,0),n_nd(ii,jj,1))=Cnds(ii,jj);
                    
               else
                    G_adj(n_nd(ii,jj,0),n_nd(ii,jj-1,0))=Gwl;
                    G_adj(n_nd(ii,jj,0),n_nd(ii,jj+1,0))=Gwl;
                    G_adj(n_nd(ii,jj,0),n_nd(ii,jj,1))=Cnds(ii,jj);
                    
               end
               
               if(ii==1)
                    G_adj(n_nd(ii,jj,1),n_nd(ii+1,jj,1))=Gbl;
                    G_adj(n_nd(ii,jj,1),n_nd(ii,jj,0))=Cnds(ii,jj);
                    
               elseif(ii==m)
                    G_adj(n_nd(ii,jj,1),n_nd(ii-1,jj,1))=Gbl;
                    G_adj(n_nd(ii,jj,1),2*m*n+m+1)=Gbl;
                    G_adj(n_nd(ii,jj,1),n_nd(ii,jj,0))=Cnds(ii,jj);
                    
               else
                    G_adj(n_nd(ii,jj,1),n_nd(ii-1,jj,1))=Gbl;
                    G_adj(n_nd(ii,jj,1),n_nd(ii+1,jj,1))=Gbl;
                    G_adj(n_nd(ii,jj,1),n_nd(ii,jj,0))=Cnds(ii,jj);
                    
               end
               
          end
          
     end
     
     for ii=1:m
          G_adj(2*m*n+ii,n_nd(ii,1,0))=Gwl;
          
     end
     for jj=1:n
          G_adj(2*m*n+m+1,n_nd(m,jj,1))=Gbl;
          
     end
     clear n_nd;

     

end


function idx=n_nd(ii,jj,fl,varargin)
     persistent n;
     if(~isempty(varargin))
          n=varargin{1};
     end
     idx=2*n*(ii-1)+2*(jj-1)+1+fl;
     

end

function [G_adj, Vin, Cnds]=init_cb(m,n,batch_size,Gwl,Gbl,fl_sym)
     if(fl_sym==1)
          Gwl=sym("Gwl");
          Gbl=sym("Gbl");
          Vin=sym("V",[m,batch_size]);
          Cnds=sym("G",[m,n]);
          
     else
          Vin=rand(m,batch_size);
          %Cnds=zeros(m,n)+10;
          Cnds=1./((10+rem(rand(m,n)*1000,991))*10e0);
          
     end
     G_adj=adj_m_w(Cnds,Gwl,Gbl); 

end

function [pruned_nds,nds_tgt]=pune_ntd(nds_td,m,n)
     n_nd(0,0,0,n);
     nds_tgt=zeros(1,n);
     for jj=1:n
          nds_tgt(jj)=n_nd(m,jj,1);
          
          
     end
     nds_td(nds_tgt)=0;
     pruned_nds=nds_td(nds_td>0);

     clear n_nd;

end

function cnds=get_mesh_cnds(cnds_inp)
     cnds=zeros((length(cnds_inp)*(length(cnds_inp)-1))/2,1);
     cnt=1;
     tot_cnds=sum(cnds_inp);
     for ii=1:length(cnds_inp)-1
          for jj=ii+1:length(cnds_inp)
               cnds(cnt)=1./((1./cnds_inp(ii))*(1./cnds_inp(jj))*tot_cnds);
               cnt=cnt+1;
               
               
               
               
          end
          
     end
     
     

end

function [row0,col0,cnds]=get_mesh_cnds_(rows,cnds_inp)
     cnds=zeros((length(cnds_inp)*(length(cnds_inp)-1)),1);
     row0=zeros((length(cnds_inp)*(length(cnds_inp)-1)),1);
     col0=zeros((length(cnds_inp)*(length(cnds_inp)-1)),1);
     cnt=1;
     tot_cnds=sum(cnds_inp);
     for ii=1:length(cnds_inp)-1
          for jj=ii+1:length(cnds_inp)
               row0(cnt)=rows(ii);
               col0(cnt)=rows(jj);
               cnds(cnt)=1./((1./cnds_inp(ii))*(1./cnds_inp(jj))*tot_cnds);
               cnt=cnt+1;
               
          end
          
     end
     row0(cnt:length(row0))=col0(1:cnt-1);
     col0(cnt:length(col0))=row0(1:cnt-1);
     cnds(cnt:length(cnds))=cnds(1:cnt-1);
     
     

end

function G1_adj=star_mesh_one_iter(G_adj,nds_td)
     nd_nums=nds_td;
     G1_adj=G_adj;
     for ii=1:length(nd_nums)
          [~,col,v]=find(G1_adj(nd_nums(ii),:));
          cnds=get_mesh_cnds(v);
          cnt=1;
          for jj=1:length(col)-1
               for kk=jj+1:length(col)
                    G1_adj(col(jj),col(kk))=G1_adj(col(jj),col(kk))+cnds(cnt);
                    G1_adj(col(kk),col(jj))=G1_adj(col(jj),col(kk));
                    cnt=cnt+1;
                    
                         
                         
                    
                    
               end
               
          end
          G1_adj(nd_nums(ii),:)=zeros(1,size(G1_adj,2));
          G1_adj(:,nd_nums(ii))=zeros(size(G1_adj,1),1);
          
          
     end
     

end

function [G_m,Ivec]=adj_to_lapl(G_adj,m,n,Vin)
     G_m=sparse(2*m*n,2*m*n);
     batch_size=size(Vin,2);
     Ivec=sparse(2*m*n,batch_size);
     for ii=1:2*m*n
          [~,col,val]=find(G_adj(ii,1:2*m*n));
          G_m(ii,ii)=sum(val)+sum(G_adj(ii,2*m*n+1:size(G_adj,2)));
          
          for jj=1:length(col)
               G_m(ii,col(jj))=-val(jj);
               
          end
          [~,col,val]=find(G_adj(ii,2*m*n+1:2*m*n+m));
          for kk=1:batch_size
               curr_acc=0;
               for jj=1:length(col)
                    curr_acc=curr_acc+val(jj)*Vin(col(jj),kk);



               end
               Ivec(ii,kk)=curr_acc;
               
          end
          
          
          
     end

end


function [Lm,Ivec]=gen_lapl(Cnds,Gwl,Gbl,Vin)
     m=size(Cnds,1); n=size(Cnds,2); batch_size=size(Vin,2);
     Lm=sparse(1,1,0,2*m*n,2*m*n,m*(6+4*(n-2))+n*(6+4*(m-2)));
     Ivec=sparse(1,1,0,2*m*n,batch_size,m*batch_size);
     n_nd(0,0,0,n);
     for ii=1:m
          for jj=1:n
               if(jj==1)
                    Lm(n_nd(ii,jj,0),n_nd(ii,jj,0))=2*Gwl+Cnds(ii,jj);
                    Lm(n_nd(ii,jj,0),n_nd(ii,jj+1,0))=-Gwl;
                    Lm(n_nd(ii,jj,0),n_nd(ii,jj,1))=-Cnds(ii,jj);
                    
               elseif(jj==n)
                    Lm(n_nd(ii,jj,0),n_nd(ii,jj,0))=Gwl+Cnds(ii,jj);
                    Lm(n_nd(ii,jj,0),n_nd(ii,jj-1,0))=-Gwl;
                    Lm(n_nd(ii,jj,0),n_nd(ii,jj,1))=-Cnds(ii,jj);
                    
               else
                    Lm(n_nd(ii,jj,0),n_nd(ii,jj,0))=2*Gwl+Cnds(ii,jj);
                    Lm(n_nd(ii,jj,0),n_nd(ii,jj-1,0))=-Gwl;
                    Lm(n_nd(ii,jj,0),n_nd(ii,jj+1,0))=-Gwl;
                    Lm(n_nd(ii,jj,0),n_nd(ii,jj,1))=-Cnds(ii,jj);
                    
               end
               
               if(ii==1)
                    Lm(n_nd(ii,jj,1),n_nd(ii,jj,1))=Gbl+Cnds(ii,jj);
                    Lm(n_nd(ii,jj,1),n_nd(ii,jj,0))=-Cnds(ii,jj);
                    Lm(n_nd(ii,jj,1),n_nd(ii+1,jj,1))=-Gbl;
                    
               elseif(ii==m)
                    Lm(n_nd(ii,jj,1),n_nd(ii,jj,1))=2*Gbl+Cnds(ii,jj);
                    Lm(n_nd(ii,jj,1),n_nd(ii,jj,0))=-Cnds(ii,jj);
                    Lm(n_nd(ii,jj,1),n_nd(ii-1,jj,1))=-Gbl;
                    
               else
                    Lm(n_nd(ii,jj,1),n_nd(ii,jj,1))=2*Gbl+Cnds(ii,jj);
                    Lm(n_nd(ii,jj,1),n_nd(ii,jj,0))=-Cnds(ii,jj);
                    Lm(n_nd(ii,jj,1),n_nd(ii-1,jj,1))=-Gbl;
                    Lm(n_nd(ii,jj,1),n_nd(ii+1,jj,1))=-Gbl;
                    
               end
               
          end
     end
     
     for ii=1:m
          for kk=1:batch_size
               Ivec(n_nd(ii,1,0),kk)=Vin(ii,kk)*Gwl;
               
          end
          
     end
     clear n_nd;

end
