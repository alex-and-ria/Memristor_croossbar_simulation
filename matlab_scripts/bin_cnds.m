
































n_dt=40;
rmse=zeros(1,n_dt);
powers=2:8;
batch_size=1;Gwl=1/10; Gbl=1/10;
R_off=1e6; R_on=1e4; V_1=0.15; V_0=0;
Gwl_swp=0.001:0.01:0.1;
Gwl_swp=1./Gwl_swp;
%out_bin=zeros(n,batch_size,length(Gwl_swp)); out_ref=zeros(n,batch_size,length(Gwl_swp)); out_wb=zeros(n,batch_size,length(Gwl_swp));
rmse_arr=zeros(length(Gwl_swp),length(powers));

for m=2.^powers
%m=5;
n=m;
diff_arr=zeros(max(n_dt+n-1,n),1);
cnt_rmse=1;
%
Vin=(rand(m,batch_size)>0.5);%rand(m,batch_size);
Cnds=(rand(m,n)>0.5);%zeros(m,n)+10;
Cnds1=(Cnds==1)*(1/R_on)+(Cnds==0)*(1/R_off);
Vin1=(Vin==1)*V_1+(Vin==0)*V_0;

%
cont_fl=0;
%for ii=1:length(Gwl_swp)
ii=1;
while ii<=length(Gwl_swp)
     Gwl=Gwl_swp(ii);
     Gbl=Gwl_swp(ii);
     G_adj=adj_m_w(Cnds,Gwl,Gbl,R_on,R_off); 
     if(cont_fl==1)
          
          Vin=(rand(m,batch_size)>0.5);
          Vin1=(Vin==1)*V_1+(Vin==0)*V_0;
          
     end

units_mn=Units(m,n,1,V_0,V_1,R_on,R_off,Gwl,Gbl);

%[G_adj, Vin, Cnds]=init_cb(m,n,batch_size,Gwl,Gbl,0,R_on,R_off);

%Cnds1=(Cnds==1)*(1/R_on)+(Cnds==0)*(1/R_off);
%Vin1=(Vin==1)*V_1+(Vin==0)*V_0;

[G_m,Ivec]=adj_to_lapl(G_adj,m,n,Vin1,batch_size);
nd_v=G_m\Ivec;

n_nd(0,0,0,n);
out_I_Gwb=zeros(n,batch_size);
for jj=1:n
     for kk=1:batch_size
          out_I_Gwb(jj,kk)=nd_v(n_nd(m,jj,1),kk)*G_adj(2*m*n+m+1,n_nd(m,jj,1));
          
     end
end
clear n_nd;

%out_I1=Cnds1'*Vin1; onr_I1=V_1/R_on;
out_I2=Cnds'*Vin;

%[out_I2,out_I1/onr_I1,out_I_Gwb./(mean(units_mn)')]

%
%out_bin(:,:,ii)=out_I2;
%out_ref(:,:,ii)=out_I1/onr_I1;
normlzr=(sum(units_mn.*Cnds))./sum(Cnds);
normlzr(isnan(normlzr))=0;
%out_wb(:,:,ii)=out_I_Gwb./normlzr';

if(n<n_dt && cnt_rmse<n_dt)
     diff_arr(cnt_rmse:cnt_rmse+n-1)=out_I_Gwb./normlzr'-out_I2;
     cnt_rmse=cnt_rmse+n;
     cont_fl=1;
     %ii=ii-1;
end
     
if(n<n_dt && cnt_rmse>=n_dt)
     diff_arr(isnan(diff_arr))=0;
     rmse_arr(ii,log(m)/log(2)-1)=sqrt(sum(diff_arr(1:cnt_rmse-1).^2)/(cnt_rmse-1));
     cnt_rmse=1;
     ii=ii+1;
     cont_fl=0;
elseif(n>=n_dt)
     rmse_arr(ii,log(m)/log(2)-1)=sqrt(sum((out_I_Gwb./normlzr'-out_I2).^2)/n);
     ii=ii+1;
end
%%figure(1); plot(1:n,out_bin(:,1:ii),'r',1:n,out_ref(:,1:ii),'g',1:n,out_wb(:,1:ii),'b');

%pause
%
end
figure(log(m)/log(2)-1); plot(1./Gwl_swp,rmse_arr(:,log(m)/log(2)-1)); title(num2str(m));
end
%for ii=1:length(powers)
%     figure(ii); plot(1./Gwl_swp,rmse_arr(:,ii)); title(num2str(2^powers(ii)));
     
     
     
     
     
%end
%figure(1); plot(Gwl_swp,out_bin,'r',Gwl_swp,out_ref,'-',Gwl_swp,out_wb,'.');







function G_adj=adj_m_w(Cnds,Gwl,Gbl,R_on,R_off)
     %node numberring is as n_nd puts them; then nodes with sources follow,
     %then ground;
     m=size(Cnds,1); n=size(Cnds,2);
     G_adj=sparse(1,1,0,2*m*n+m+1,2*m*n+m+1,m*(5+3*(n-2))+n*(5+3*(m-2))+m+1);
     n_nd(0,0,0,n);
     Cnds1=(Cnds==1)*(1/R_on)+(Cnds==0)*(1/R_off);
     for ii=1:m
          for jj=1:n
               if(jj==1)
                    G_adj(n_nd(ii,jj,0),2*m*n+ii)=Gwl;
                    G_adj(n_nd(ii,jj,0),n_nd(ii,jj+1,0))=Gwl;
                    G_adj(n_nd(ii,jj,0),n_nd(ii,jj,1))=Cnds1(ii,jj);
                    
               elseif(jj==n)
                    G_adj(n_nd(ii,jj,0),n_nd(ii,jj-1,0))=Gwl;
                    G_adj(n_nd(ii,jj,0),n_nd(ii,jj,1))=Cnds1(ii,jj);
                    
               else
                    G_adj(n_nd(ii,jj,0),n_nd(ii,jj-1,0))=Gwl;
                    G_adj(n_nd(ii,jj,0),n_nd(ii,jj+1,0))=Gwl;
                    G_adj(n_nd(ii,jj,0),n_nd(ii,jj,1))=Cnds1(ii,jj);
                    
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

function [G_adj, Vin, Cnds]=init_cb(m,n,batch_size,Gwl,Gbl,fl_sym,R_on,R_off)
     if(fl_sym==1)
          Gwl=sym("Gwl");
          Gbl=sym("Gbl");
          Vin=sym("V",[m,batch_size]);
          Cnds=sym("G",[m,n]);
          
     else
          %Vin=ones(m,1);%(rand(m,batch_size)>0.5);%rand(m,batch_size);
          %Cnds=ones(m,n);%(rand(m,n)>0.5);%zeros(m,n)+10;\
          Vin=(rand(m,batch_size)>0.5);%rand(m,batch_size);
          Cnds=(rand(m,n)>0.5);%zeros(m,n)+10;
          
          %Cnds=1./((10+rem(rand(m,n)*1000,991))*10e0);
          
     end
     G_adj=adj_m_w(Cnds,Gwl,Gbl,R_on,R_off); 

end

function [G_m,Ivec]=adj_to_lapl(G_adj,m,n,Vin,batch_size)
     G_m=sparse(2*m*n,2*m*n);
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


function units_mn=Units(m,n,fl_wb,V_0,V_1,R_on,R_off,Gwl,Gbl)
     Cnds=zeros(m,n);
     Vin=zeros(m,1);
     units_mn=zeros(m,n);
     for ii=1:m
          Vin(ii,1)=1;
          Vin1=(Vin==1)*V_1+(Vin==0)*V_0;
          for jj=1:n
               Cnds(ii,jj)=1;
               
               if(fl_wb==1)
                    G_adj=adj_m_w(Cnds,Gwl,Gbl,R_on,R_off);
                    [G_m,Ivec]=adj_to_lapl(G_adj,m,n,Vin1,1);
                    nd_v=G_m\Ivec;

                    n_nd(0,0,0,n);
                    units_mn(ii,jj)=nd_v(n_nd(m,jj,1))*G_adj(2*m*n+m+1,n_nd(m,jj,1));
                    
                    clear n_nd;
                    
               else
                    Cnds1=(Cnds==1)*(1/R_on)+(Cnds==0)*(1/R_off);
                    out_I=Cnds1'*Vin1;
                    units_mn(ii,jj)=out_I(jj,1);
                    
                    
               end
               Cnds(ii,jj)=0;
               
          end
          Vin(ii,1)=0;

     end
    
end
