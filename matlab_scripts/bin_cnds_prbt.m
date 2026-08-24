
































pwr0=3;powers=pwr0:5;
batch_size=1;Gwl=1/10; Gbl=1/10;
R_off=1e6; R_on=1e4; V_1=0.15; V_0=0;
Gwl_swp=1e-5:5:20;
Gwl_swp=1./Gwl_swp;




% for m=4%2.^powers
%      m=4;
%      n=m;
% %      m=2;
% %      n=2;
% %      V_1=5; V_0=5;
% %      R_on=1000; R_off=1000;
%      Cnds=(rand(m,n)>0.5);%zeros(m,n)+10;
%      Cnds1=(Cnds==1)*(1/R_on)+(Cnds==0)*(1/R_off);
%      Vin=(rand(m,batch_size)>0.5);%rand(m,batch_size);
%      Vin1=(Vin==1)*V_1+(Vin==0)*V_0;
% 
%      out_I_Gwb=zeros(n,batch_size,length(Gwl_swp));
%      out_wb=zeros(n,batch_size,length(Gwl_swp));
% 
%      for(ii=1:length(Gwl_swp))
%           Gwl=Gwl_swp(ii);
%           Gbl=Gwl_swp(ii);
%           G_adj=adj_m_w(Cnds,Gwl,Gbl,R_on,R_off);
%           Vin=(rand(m,batch_size)>0.5);
%           Vin1=(Vin==1)*V_1+(Vin==0)*V_0;
% 
%           units_mn=Units(m,n,1,V_0,V_1,R_on,R_off,Gwl,Gbl);
%           [G_m,Ivec]=adj_to_lapl(G_adj,m,n,Vin1,batch_size);
%           nd_v=G_m\Ivec;
% 
%           n_nd(0,0,0,n);
%           for jj=1:n
%                for kk=1:batch_size
%                     out_I_Gwb(jj,kk,ii)=nd_v(n_nd(m,jj,1),kk)*G_adj(2*m*n+m+1,n_nd(m,jj,1));
% 
%                end
%           end
%           clear n_nd;
% 
%           out_wb_ref=Cnds'*Vin;
%           I1_ref=V_1/R_on;
%           out_wb1=Cnds1'*Vin1/I1_ref;
% 
%           normlzr=(sum(units_mn.*Cnds))./sum(Cnds);
%           normlzr(isnan(normlzr))=0;
%           out_wb(:,:,ii)=out_I_Gwb(:,:,ii)./normlzr';
%           figure(1); plot(1:n,out_wb_ref,'r',1:n,out_wb1,'g',1:n,out_wb(:,:,ii),'b',1:n,out_I_Gwb(:,:,ii)./I1_ref,'b--');
% 
%      end
% end


ii0_min=zeros(length(powers),1);
jj_rec=zeros(length(powers),1);
ii_rec=zeros(length(powers),1);
for m=2.^powers
     n=m;
     Wl_Bl_min=0.01; Wl_Bl_maz=100; acc=1-2;
     Vin1=ones(m,1)*V_1;
%           Cmbn=zeros(2^m-1,m);
%           for(ii=1:2^m-1)
%                Cmbn(ii,:)=dec2bin(ii,m)-'0';
% 
%           end
%           sort_sums=sum(Cmbn,2);
%           [sort_sums,idx]=sort(sort_sums);
%           Cmbn=Cmbn(idx,:);
%           indxs=zeros(m,1);%indxs=zeros(m+1,1);
%           indxs(1)=1;
%      %      idx_cnt=1;
%      %      for(ii=2:2^m-1)
%      %           if(sort_sums(ii)~=sort_sums(ii-1))
%      %                idx_cnt=idx_cnt+1;
%      %                indxs(idx_cnt)=ii;
%      %           end
%      %           
%      %      end%effectivelly cummulative sum in increments of nchoosek(m,ii);
% 
%           for(ii=2:m)
%                indxs(ii)=indxs(ii-1)+nchoosek(m,ii-1);
% 
%           end
%           %indxs(m+1)=2^m;
%           %currents=cell(n,m,length(Gwl_swp));%{output,number,line_cnds}
%      
     currents_max=zeros(n,m,length(Gwl_swp));
     currents_min=zeros(n,m,length(Gwl_swp));
     for(jj=1:n)%output column;
          %Cnds=zeros(m,n);
          for(ii=1:m)%number to output in column;
               %curr_curents=zeros(nchoosek(m,ii),1);
               %for minimum current in j_th ouput line, of number ii;
               Cnds=ones(m,n);
               Cnds(:,jj)=zeros(m,1);
               for(ii0=1:ii)
                    Cnds(ii0,jj)=1;
               end
               for(ii0=1:length(Gwl_swp))
                    Gwl=Gwl_swp(ii0);
                    Gbl=Gwl_swp(ii0);
                    G_adj=adj_m_w(Cnds,Gwl,Gbl,R_on,R_off);
                    [G_m,Ivec]=adj_to_lapl(G_adj,m,n,Vin1,batch_size);
                    nd_v=G_m\Ivec;

                    n_nd(0,0,0,n);
                    currents_min(jj,ii,ii0)=nd_v(n_nd(m,jj,1))*G_adj(2*m*n+m+1,n_nd(m,jj,1));
                    clear n_nd;
                    
               end
               %for maximum current in j_th ouput line, of number ii;
               Cnds=zeros(m,n);
               Cnds(:,jj)=zeros(m,1);
               for(ii0=1:ii)
                    Cnds(m-ii0+1,jj)=1;
               end
               for(ii0=1:length(Gwl_swp))
                    Gwl=Gwl_swp(ii0);
                    Gbl=Gwl_swp(ii0);
                    G_adj=adj_m_w(Cnds,Gwl,Gbl,R_on,R_off);
                    [G_m,Ivec]=adj_to_lapl(G_adj,m,n,Vin1,batch_size);
                    nd_v=G_m\Ivec;

                    n_nd(0,0,0,n);
                    currents_max(jj,ii,ii0)=nd_v(n_nd(m,jj,1))*G_adj(2*m*n+m+1,n_nd(m,jj,1));
                    clear n_nd;
                    
               end
                    
                    
                    
%                     Cnds=zeros(m,n);
%                     %Cnds=ones(m,n);
%                     %Cnds(:,jj)=zeros(m,1);
%                     curr_curents=zeros(nchoosek(m,ii),length(Gwl_swp));
%                     for(ii0=1:length(Gwl_swp))
%                          for(kk=1:nchoosek(m,ii))
%                               Cnds(:,jj)=Cmbn(indxs(ii)+kk-1,:);
%                               Gwl=Gwl_swp(ii0);
%                               Gbl=Gwl_swp(ii0);
%                               G_adj=adj_m_w(Cnds,Gwl,Gbl,R_on,R_off);
%                               [G_m,Ivec]=adj_to_lapl(G_adj,m,n,Vin1,batch_size);
%                               nd_v=G_m\Ivec;
%                               
%                               n_nd(0,0,0,n);
%                               curr_curents(kk,ii0)=nd_v(n_nd(m,jj,1))*G_adj(2*m*n+m+1,n_nd(m,jj,1));
%                               clear n_nd;
%                               
%                          end
% 
%                     end
               
               
               
          end
          
          
     end
     iter=log2(m)-pwr0+1;
     ii0_min(iter)=length(Gwl_swp);
     jj_rec(iter)=0;
     ii_rec(iter)=0;
     ovlp_fl=0;
     for(jj=1:n)
          for(ii=1:m-1)
               for(ii0=1:ii0_min(iter))
                    if(currents_max(jj,ii,ii0)>currents_min(jj,ii+1,ii0))
                         ii0_min(iter)=ii0;
                         jj_rec(iter)=jj;
                         ii_rec(iter)=ii;
                         ovlp_fl=1;
                         break;
                         
                    end
                    
               end
               
          end
          
     end
     if(ovlp_fl==1)
          ii0_min(iter)=ii0_min(iter)-1;%recoreded smallest indexwith overlap, so 
          %now need to decrement it, to have nonoverlaping regin;
     end
     
     
     
     close all
     for(jj=1:n)
          figure(jj);
          %x_ax_data=repmat(1:m,1,length(Gwl_swp));
          %plot(x_ax_data,currents_min(jj,:),'.',x_ax_data,currents_max(jj,:),'*',x_ax_data,V_1/R_on*x_ax_data,'+')
          x_ax_data=zeros(2,length(Gwl_swp)*m);
          for(kk=1:m)
               if(kk==1)
                    idx_l=1;
               else
                    idx_l=(kk-1)*length(Gwl_swp)+1;
                    
               end
               x_ax_data(:,idx_l)=kk;
               for(kk0=1:length(Gwl_swp)-1)%when kk0==0, it is actually the first set of currents (min and max);
                    x_ax_data(:,idx_l+kk0)=...
                         x_ax_data(:,idx_l+kk0-1)+1/(length(Gwl_swp)+1);%shifts for better visualization;
                    
                    
               end
               
          end
          
          
          curr_min=squeeze(currents_min(jj,:,:))';
          curr_max=squeeze(currents_max(jj,:,:))';
          currents_rng=[curr_min(:)'; curr_max(:)'];
          
          hold on
          plot(x_ax_data,currents_rng,'.-');
          
          x_ax_ref=zeros(1,2*m);
          y_ax_ref=zeros(1,2*m);
          for(kk=0:m-1)
               x_ax_ref(2*kk+1)=kk+0.5;
               x_ax_ref(2*kk+2)=kk+1.5;
               y_ax_ref(2*kk+1:2*kk+2)=kk+1;
               
          end
          
          plot(x_ax_ref,y_ax_ref*V_1/R_on);
          hold off
          
          %plot(x_ax_data,currents_rng,'.-',);'*',x_ax_data,V_1/R_on*repmat(1:m,length(Gwl_swp),1),'+')
          
          
          if(ii0_min(iter)>0)
               figure(n+jj);
               curr_min=squeeze(currents_min(jj,:,ii0_min(iter)))';
               curr_max=squeeze(currents_max(jj,:,ii0_min(iter)))';
               currents_rng=[curr_min(:)'; curr_max(:)'];

               hold on
               plot(repmat(1:m,2,1),currents_rng);
               plot(x_ax_ref,y_ax_ref*V_1/R_on);
               hold off
               
               
          end
          

          
     end
     
     
     
     
     
end

figure(1); plot(1:legth(powers),1/Gwl_swp(ii0_min));

firgure(2); 



% for m=2.^powers
%      n=m;
%      Vin1=ones(m,1)*V_1;
% %      for jj=n-1:-1:0
% %           range=range+2^jj*m;
% %      end%effectively m*(2^n-1);
%      range=m*(2^n-1);
%      min_crnt=zeros(range,n);
%      max_crnt=zeros(range,n);
%      for(out_val=1:range)
%           Cnds=zeros(m,n);
%           val_curr_min=out_val;
%           jj=n;
%           while(val_curr_min>0)%TODO redo 
%                rem_curr=rem(val_curr_min,m+1);
%                kk=1;
%                while(rem_curr>0)
%                     Cnds(kk,jj)=1;
%                     rem_curr=rem_curr-1; kk=kk+1;
%                     
%                end
%                jj=jj-1;
%                val_curr_min=fix(val_curr_min/4);
%                
%                
%               
%                
%           end
%           G_adj=adj_m_w(Cnds,Gwl,Gbl,R_on,R_off);
%           [G_m,Ivec]=adj_to_lapl(G_adj,m,n,Vin1,batch_size);
%           nd_v=G_m\Ivec;
% 
%           n_nd(0,0,0,n);
%           for jj=1:n
%                min_crnt(out_val,jj)=nd_v(n_nd(m,jj,1))*G_adj(2*m*n+m+1,n_nd(m,jj,1));
%           
%           end
%           clear n_nd;
%           
%           
%           
%           Cnds=zeros(m,n);
%           val_curr_max=out_val;
%           prev_jj=n;
%           while(val_curr_max>0)
%                jj=fix(log2(val_curr_max));
%                if(jj>(prev_jj-1))
%                     jj=prev_jj-1;
%                end
%                rem_curr=fix(val_curr_max/2^jj);
%                if(rem_curr>m)
%                     rem_curr=m;
%                end
%                base_curr=val_curr_max-rem_curr*2^jj;
%                while(rem_curr>0)
%                     kk=m;
%                     while(rem_curr>0)
%                          Cnds(kk,n-jj)=1;
%                          rem_curr=rem_curr-1; kk=kk-1;
% 
%                     end
% 
% 
% 
%                end
%                val_curr_max=base_curr;
%                prev_jj=prev_jj-1;
%                
%           end
%           
%           
%           
%           
%      end
%      
% end

% 
% for m=2.^powers
%      n=m;
%      Cnds1=ones(m,n)*1/R_on;
%      Cmbn=zeros(2^m-1,m);
%      for(ii=1:2^m-1)
%           Cmbn(ii,:)=dec2bin(ii,m)-'0';
%           
%      end
%      sort_sums=sum(Cmbn,2);
%      [sort_sums,idx]=sort(sort_sums);
%      Cmbn=Cmbn(idx,:);
%      indxs=zeros(m+1,1);
%      indxs(1)=1;
% %      idx_cnt=1;
% %      for(ii=2:2^m-1)
% %           if(sort_sums(ii)~=sort_sums(ii-1))
% %                idx_cnt=idx_cnt+1;
% %                indxs(idx_cnt)=ii;
% %           end
% %           
% %      end%effectivelly cummulative sum in increments of nchoosek(m,ii);
%      for(ii=2:m)
%           indxs(ii)=indxs(ii-1)+nchoosek(m,ii-1);
%           
%      end
%      %idx_cnt=idx_cnt+1;
%      indxs(m+1)=2^m;
% %    for jj=n-1:-1:0
% %         range=range+2^jj*m;
% %    end%effectively m*(2^n-1);
%      range=m*(2^n-1);
%      tgt_val=zeros(n,1);
%      tgt_combs=cell(n,1);
%      for(out_val=1:range)
%           val0=0;
%           rems=out_val;
%           jj=n; cnt_pw=0;
%           %if(val0 && rems are correct)
%           while(rems>0)
%                while(rems>0)
%                     tgt_val(jj)=rem(rems,m+1);
%                     rems=fix(rems/m);
%                     if(tgt_val(jj)~=0)
%                          tgt_combs{jj}=Cmbn(indxs(tgt_val(jj)):indxs(tgt_val(jj)+1)-1,:);
% 
%                     end
%                     cnt_pw=cnt_pw+1;
%                     jj=n-cnt_pw;
%                     %tgt_cmbs;
% 
% 
%                end
%                while(val0>0)
%                     tgt_val(jj)=rem(val0,m+1);
%                     val0=fix(val0/m);
%                     if(tgt_val(jj)~=0)
%                          tgt_combs{jj}=Cmbn(indxs(tgt_val(jj)):indxs(tgt_val(jj)+1)-1,:);
% 
%                     end
%                     cnt_pw=cnt_pw+1;
%                     jj=n-cnt_pw;
% 
% 
% 
%                end
%                val0=val0+1;
%                rems=rem(out_val/2^val0);
%                
%           end
%          
%           
%           val0=
%           val1=fix(out_val/2^kk);
%           rems=rem(out_val/2^kk);
%           
%           while(curr_val0>0)
%                curr_val=curr_val0;
%                for(jj=1:n)
%                     tgt_val(jj)=rem(curr_val,m+1);
%                     curr_val=fix(curr_val/m);
%                     if(tgt_val(jj)~=0)
%                          tgt_combs{jj}=Cmbn(indxs(tgt_val(jj)):indxs(tgt_val(jj)+1)-1,:);
% 
%                     end
%                     %tgt_combs{jj}=Cmbn(indxs(tgt_val(jj)
%                     %for(kk=1:
% 
% 
%                end
%                rem0=(curr_val0,2);
%                curr_val0=fix(curr_val0/2);
%                
% 
% 
% 
%           end
%           
%      end
%      
%      for(ii=1:length(Gwl_swp))
%           Gwl=Gwl_swp(ii);
%           Gbl=Gwl_swp(ii);
%           G_adj=adj_m_w(Cnds,Gwl,Gbl,R_on,R_off);
%           
% %           range=0;
% %           for jj=n-1:-1:0
% %                range=range+2^jj*m;
% %           end%effectively m*(2^n-1);
%           range=m*(2^n-1);
%           tgt_val=zeros(n,1);
%           
%           for(out_val=1:range)
%                curr_val=out_val;
%                for(jj=1:n)
%                     tgt_val(jj)=rem(curr_val,m);
%                     curr_val=fix(curr_val/m);
%                     
%                     
%                end
%                comb=nchoosek(ones(n,1),tgt_val(jj));
%                
%                
%                
%           end
%           
%           
%      end
%      
%      
%      
% end





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
                    G_adj(n_nd(ii,jj,1),n_nd(ii,jj,0))=Cnds1(ii,jj);
                    
               elseif(ii==m)
                    G_adj(n_nd(ii,jj,1),n_nd(ii-1,jj,1))=Gbl;
                    G_adj(n_nd(ii,jj,1),2*m*n+m+1)=Gbl;
                    G_adj(n_nd(ii,jj,1),n_nd(ii,jj,0))=Cnds1(ii,jj);
                    
               else
                    G_adj(n_nd(ii,jj,1),n_nd(ii-1,jj,1))=Gbl;
                    G_adj(n_nd(ii,jj,1),n_nd(ii+1,jj,1))=Gbl;
                    G_adj(n_nd(ii,jj,1),n_nd(ii,jj,0))=Cnds1(ii,jj);
                    
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
