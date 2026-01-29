



























m=256; n=256;
batch_size = 1;
%Vin=rand(m,batch_size);
Vin=(zeros(m,batch_size)+5);
Gwl=1./100;
Gbl=4./100;
% conductances
%approximately 10K Ohm to 1000k Ohm;
Cnds=1./((10+(rem(rand(m,n)*1000,991)))*10e3);
%Cnds=1./((10+rem(rand(m,n)*1000,991))*10e0);
%Cnds=1./(zeros(3,3)+1000);
[Lm,Ivec]=gen_lapl(Cnds,Gwl,Gbl,Vin);
tic(); [L,U,P] = lu(Lm); toc();
tic(); y = L\(P*Ivec); toc();
tic(); x = U\y; toc();


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