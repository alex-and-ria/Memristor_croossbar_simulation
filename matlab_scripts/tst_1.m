




















fl_nm=strcat('/share/share0/abc/star_mesh/rw_.h');
[f_id,err_msg]=fopen(fl_nm,'w');
if(f_id<0)
     disp('f_id<0');
end
m=4; n=4; batch_size=1;
Gwl=1./100; Gbl=4./100;
[G_adj, Vin]=init_cb(m,n,batch_size,Gwl,Gbl,0);
[row,col,val]=find(G_adj);

nds_td=1:2*m*n;
[nds_td,nds_tgt]=pune_ntd(nds_td,m,n);

fprintf(f_id,"unsigned int row[]={%u",row(1));
for ii=2:size(row,1)
     fprintf(f_id,',%u',row(ii));

     
end
fprintf(f_id,"};\n");

fprintf(f_id,"unsigned int col[]={%u",col(1));
for ii=2:size(col,1)
     fprintf(f_id,',%u',col(ii));

     
end
fprintf(f_id,"};\n");

fprintf(f_id,"double val[]={%.15g",val(1));
for ii=2:size(val,1)
     fprintf(f_id,',%.15g',val(ii));

     
end
fprintf(f_id,"};\n");

fprintf(f_id,"unsigned int nds_td[]={%u",nds_td(1));
for ii=2:size(nds_td,2)
     fprintf(f_id,',%u',nds_td(ii));

     
end
fprintf(f_id,"};\n");

fprintf(f_id,"unsigned int len1=%u; unsigned int nds_n=%d;\n",size(row,1),size(nds_td,2));

fprintf(f_id,"unsigned int nds_td1[]={%u",nds_tgt(1));
for ii=2:size(nds_tgt,2)
     fprintf(f_id,',%u',nds_tgt(ii));

     
end
fprintf(f_id,"};\n");

fprintf(f_id,"unsigned int nds_n1=%u; unsigned int max_mx_sz=2;\n", size(nds_tgt,2));

fclose(f_id);



%calllib('libnode_schr','dense_rdct',row,rw,col,cl,val,vl,lenp,ln,nds_td_p,nds_np);


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
          Cnds=zeros(m,n)+10;
          %Cnds=1./((10+rem(rand(m,n)*1000,991))*10e0);
          
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