

































loadlibrary('libnode_schr','node_schr.h')
libfunctions('libnode_schr','-full')
batch_size=1;
fl_nm=strcat('/share/share0/abc/star_mesh/dts_one_iter_openmp',num2str(batch_size));
[f_id,err_msg]=fopen(fl_nm,'w');
if(f_id<0)
     disp('f_id<0');
end
for ii=1:8
     m=2^(ii); n=m;
     fprintf(f_id,'\n%ux%u\n',m,n);
     Gwl=1./100; Gbl=4./100;
     [G_adj, Vin]=init_cb(m,n,batch_size,Gwl,Gbl,0);
     [row,col,val]=find(G_adj);
     
     nds_td=1:2*m*n;
     nds_td=pune_ntd(nds_td,m,n);
     
     nds_td_p=libpointer('uint32Ptr',nds_td);
     tic();
     calllib('libnode_schr','node_analyzer',row,col,val,size(col,1),nds_td_p,size(nds_td,2));
     dt=toc(); fprintf(f_id,'%.15g, ',dt);
     nds_td=nds_td_p.Value(nds_td_p.Value>0);
     rw=libpointer('uint32PtrPtr'); cl=libpointer('uint32PtrPtr'); vl=libpointer('doublePtrPtr'); ln=100;%libpointer('uint32Ptr');
     
     tic()
     [~,~,~,~,~,~,~,ln,~,~]=calllib('libnode_schr','star_mesh_base',row,rw,col,cl,val,vl,size(col,1),ln,nds_td,size(nds_td,2));
     dt1=toc(); fprintf(f_id,'%.15g, ',dt1);
     rw.Value.setdatatype('uint32Ptr',1,ln); cl.Value.setdatatype('uint32Ptr',1,ln); vl.Value.setdatatype('doublePtr',1,ln);
     Gm_tf=sparse(rw.Value,cl.Value,vl.Value);
     
     tic()
     G_one_iter=star_mesh_one_iter(G_adj,nds_td);
     dt0=toc(); fprintf(f_id,'%.15g\n',dt0);
     sum(sum(abs(Gm_tf-Gm_tf)))
     
     
     
end
fclose(f_id);
unloadlibrary libnode_schr




function G_adj=adj_m_w(Cnds,Gwl,Gbl)
     %node numberring is as n_nd puts them; then nodes with sources follow,
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

function [G_adj, Vin]=init_cb(m,n,batch_size,Gwl,Gbl,fl_sym)
     if(fl_sym==1)
          Gwl=sym("Gwl");
          Gbl=sym("Gbl");
          Vin=sym("V",[m,batch_size]);
          Cnds=sym("G",[m,n]);
          
     else
          Vin=rand(m,batch_size);
          Cnds=1./((10+rem(rand(m,n)*1000,991))*10e0);
          
     end
     G_adj=adj_m_w(Cnds,Gwl,Gbl); 

end

function pruned_nds=pune_ntd(nds_td,m,n)
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