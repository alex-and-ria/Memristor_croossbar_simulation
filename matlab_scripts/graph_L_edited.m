































alex
m=2;n=3;
batch_size=2;
Gwl=1./100;
%Gwl=sym("Gwl");
Gbl=4./100;
%Gbl=sym("Gbl");
Vin=rand(m,batch_size);
%Vin=(zeros(m,batch_size)+10);
%Vin=sym("V",[m,batch_size]);
Cnds=1./((10+rem(rand(m,n)*1000,991))*10e0);
%Cnds=1./(zeros(m,n)+1000);
%Cnds=sym("G",[m,n]);

[Lm,Ivec]=gen_lapl(Cnds,Gwl,Gbl,Vin);
G_adj=adj_m_w(Cnds,Gwl,Gbl); spy(G_adj);
%[Lm1,Ivec1]=adj_to_lapl(G_adj,Vin);
%Lm\Ivec(:,1)-Lm1\Ivec1(:,1);
prot_nds=[1,2,3,4];
G1_adj=star_mesh(G_adj, prot_nds,m);
[Lm1,Ivec1]=adj_to_lapl(G1_adj,Vin);

x1=pinv(Lm1)*Ivec1; x0=pinv(Lm)*Ivec;
sum(abs(x1(prot_nds,:)-x0(prot_nds,:)),'all')
%diff=pinv(Lm1)*Ivec1-pinv(Lm)*Ivec;  
%sum((Lm\Ivec-Lm1\Ivec1)>10e-10,'all');
%diff(isnan(diff))=0;




%%
powers=1:8;
for batch_size=10.^[0:1]
%batch_size=1;
Gwl=1./100; Gbl=4./100;
fl_nm=strcat('/share/share0/abc/star_mesh/dts',num2str(batch_size));
[f_id,err_msg]=fopen(fl_nm,'w');
if(f_id<0)
     disp('f_id<0');
end
kft=2;
for m=64%2.^powers
     for n=64%2.^powers
          fprintf(f_id,'Crossbar: %ux%u\n',m,n);
          [G_adj, Vin]=init_cb(m,n,batch_size,Gwl,Gbl,0);
          n_nd(0,0,0,n);
          for p=32%2.^[1:(log(m)/log(2))]
               for q=p%2.^[1:(log(n)/log(2))]
                    partitions=zeros((m/p)*(n/q),kft*p*q);
                    ii=1; ii0=1; jj=1; jj0=1;
                    for kk=1:(m/p)*(n/q)
                         for ll=1:p*q
                              for qq=1:kft
                                   partitions(kk,kft*(ll-1)+qq)=n_nd(ii,jj,qq-1);
                                   
                              end
                              jj=jj+1;
                              if((jj-jj0)==q)%reached the end of the line in current block;
                                   jj=jj0;
                                   ii=ii+1;
                                   
                                   
                              end
                              %ii reaches end of block automatically with the end of the loop;
                              
                         end
                         jj0=jj0+q;
                         if(jj0==n+1)%reached the end of the block;
                              ii0=ii0+p;
                              jj0=1;
                              
                              
                         end
                         ii=ii0;
                         jj=jj0;
                         
                              
                         
                    end
                    
                    
                    
                    fprintf(f_id,'partitions: %ux%u\n',p,q);
                    %fprintf(f_id,',');
                    for kk=1:(m/p)*(n/q)
                         fprintf(f_id,'%u,',kk);
                         
                    end
                    fprintf(f_id,'\n');
                    dts=zeros((m/p)*(n/q),9);
                    for kk=1:(m/p)*(n/q)
                         tic(); G1_adj=star_mesh(G_adj, partitions(kk,:),m); dts(kk,1)=toc();
                         tic(); [Lm1,Ivec1]=adj_to_lapl(G1_adj,Vin); dts(kk,2)=toc();
                         Lm_q=Lm1(:,any(Lm1));
                         tic(); x0=Lm_q\Ivec1; dts(kk,3)=toc();
                         tic(); [L,U,P]=lu(Lm_q); dts(kk,4)=toc();
                         tic(); y=L\P*Ivec1; dts(kk,5)=toc();
                         tic(); x00=U\y; dts(kk,6)=toc();
                         
                         tic(); [L,U]=lu(Lm_q); dts(kk,7)=toc();
                         tic(); y=L\Ivec1; dts(kk,8)=toc();
                         tic(); x00=U\y; dts(kk,9)=toc();
                         
                         %tic(); x1=pinv(Lm1)*Ivec1; dts(kk,7)=toc();
                         
                         eps=10e-9;
                         %sum(x0(abs(x0)>eps)-x1(abs(x1)>eps),'all')+sum((x0-x00),'all')
                         
                         %add recording to file;
                         
                    end
                    
                    for qq=1:size(dts,2)
                         for kk=1:(m/p)*(n/q)
                              
                              fprintf(f_id,'%.15g,',dts(kk,qq));
                              
                         end
                         fprintf(f_id,'\n');
                         
                         
                    end
                    fprintf(f_id,'totals: transf. time=%.15g; total LU=%.15g; backward+forward substitutions=%.15g; pinv=%.15g\n',sum(dts(:,1)),sum(dts(:,4)),sum(dts(:,5))+sum(dts(:,6)),sum(dts(:,7)));
                    fprintf(f_id,'\n');
                    
                         
               end
             
                    
               
          end
          clear n_nd;
          fprintf(f_id,'\n');
          
     end
     
end

fclose(f_id);
end

%%


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

function [Lm,Ivec]=adj_to_lapl(G_adj,Vin)
     m=size(Vin,1); batch_size=size(Vin,2);
     %n=(size(G_adj,2)-1-m)/(2*m);
     %Lm=sparse(1,1,0,2*m*n,2*m*n,m*(6+4*(n-2))+n*(6+4*(m-2)));
     %Lm=sparse(size(G_adj,1)-1-m,size(G_adj,2)-1-m);
     %Ivec=sparse(1,1,0,2*m*n,batch_size,m*batch_size);
%      n_nd(0,0,0,n);
%      
%      for ii=1:m
%           for jj=1:n
%                [~,col,v]=find(G_adj(n_nd(ii,jj,0),1:2*m*n));
%                Lm(n_nd(ii,jj,0),n_nd(ii,jj,0))=sum(G_adj(n_nd(ii,jj,0),:));
%                for kk=1:length(col)
%                     Lm(n_nd(ii,jj,0),col(kk))=-v(kk);
%                     
%                end
%                [~,col,v]=find(G_adj(n_nd(ii,jj,1),1:2*m*n));
%                Lm(n_nd(ii,jj,1),n_nd(ii,jj,1))=sum(G_adj(n_nd(ii,jj,1),:));
%                for kk=1:length(col)
%                     Lm(n_nd(ii,jj,1),col(kk))=-v(kk);
%                     
%                end
%                
%                
%           end
%           
%      end
     
     %convection is that current goes from node with bigger number to the
     %node with smaller number; as long as currents taken into account
     %consistently, calculations should be fine;
     for ii=1:(size(G_adj,1)-1-m)
          [~,col,v]=find(G_adj(ii,1:(size(G_adj,2)-1-m)));
          Lm(ii,ii)=sum(G_adj(ii,:));
          for jj=1:length(col)
               Lm(ii,col(jj))=-v(jj);
               
          end
          
     end
     
     m1=size(G_adj,1)-m-1;%nodes that are not sources or ground;
     for ii=1:m1
          [~,col,~]=find(G_adj(ii,m1+1:(size(G_adj,1)-1)));%col vector starts from 1 which is m1+1 column for G_adj;
          for kk=1:batch_size
               curr_acc=0;
               for jj=1:length(col)
                    curr_acc=curr_acc+G_adj(ii,m1+col(jj))*Vin(col(jj),kk);%warning: can depend on numbering scheme;
                    
                    
               end
               Ivec(ii,kk)=curr_acc;
               
          end
          
     end
     
%      for ii=1:m
%           for kk=1:batch_size
%                Ivec(n_nd(ii,1,0),kk)=Vin(ii,kk)*G_adj(n_nd(ii,1,0),2*m*n+ii);
%                
%           end
%           
%      end
%      
%      clear n_nd;

end

function G1_adj=star_mesh(G_adj, prot_nds,n_srcs)
     %nd_nums=1:size(G_adj,1)-n_srcs-1;
     %tic();
     %nd_nums(ismember(nd_nums,prot_nds))=[];%delete protected nodes form list of nodes list;
     %toc();
     nd_nums=1:size(G_adj,1)-n_srcs-1;
     %tic();
     for ii=1:length(nd_nums)
          for jj=1:length(prot_nds)
               if(nd_nums(ii)==prot_nds(jj))
                    nd_nums(ii)=0;
                    
               end
               
          end
          
     end
     nd_nums(nd_nums==0)=[];
     %toc();
     
     G1_adj=G_adj;
     for ii=1:length(nd_nums)
          [~,col,v]=find(G1_adj(nd_nums(ii),:));
          cnds=get_mesh_cnds(v);
          cnt=1;
          for jj=1:length(col)-1
               for kk=jj+1:length(col)
%                     if(G1_adj(col(jj),col(kk))~=0)
%                          G1_adj(col(jj),col(kk))=G1_adj(col(jj),col(kk))+cnds(cnt);
                          %if there is already resistor, use rule of
                          %parrallel resistors connection;

%                     else
%                          G1_adj(col(jj),col(kk))=cnds(cnt);
%                     end
                    
                    G1_adj(col(jj),col(kk))=G1_adj(col(jj),col(kk))+cnds(cnt);
                    %G1_adj(col(jj),col(kk))=sym("G_eq");
                    G1_adj(col(kk),col(jj))=G1_adj(col(jj),col(kk));
                    cnt=cnt+1;
                    
                         
                         
                    
                    
               end
               
          end
          G1_adj(nd_nums(ii),:)=zeros(1,size(G1_adj,2));
          G1_adj(:,nd_nums(ii))=zeros(size(G1_adj,1),1);
          
          
     end
     

end

function cnds=get_mesh_cnds(cnds_inp)
     cnds=zeros((length(cnds_inp)*(length(cnds_inp)-1))/2);
     cnt=1;
     tot_cnds=sum(cnds_inp);
     for ii=1:length(cnds_inp)-1
          for jj=ii+1:length(cnds_inp)
               cnds(cnt)=1./((1./cnds_inp(ii))*(1./cnds_inp(jj))*tot_cnds);
               cnt=cnt+1;
               
               
               
               
          end
          
     end
     
     

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

