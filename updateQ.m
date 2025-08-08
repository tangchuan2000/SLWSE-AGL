function [Q] = updateQ(S, Q0)
% Input
% S   anchor graph n*m
% Q0  is the initial label matrix m*c


[n,m] = size(S);
Q = Q0;
%% store once
aa=sum(Q,1);
[~,label]=max(Q,[],2);
SSS=2*(S'*S);
XX=diag(SSS)./2;

SSQQ= SSS* Q;% SSQQ(i,:)   bh'Sul
ybby=diag(Q'*SSQQ/2);
%% compute Initial objective function value
  [~,Label0]=max(S*Q0,[],2);
  obj(1)= sum((ybby).^(1/2));
%%
for iter=1:30  
     for i = 1:m   
         mm = label(i) ;
        if aa(mm)==1
            continue;  
        end    
       %% 
        V2=ybby'+(SSQQ(i,:)+XX(i)).*(1-Q(i,:));
        V1=(ybby'-(SSQQ(i,:)-XX(i)).*Q(i,:));
        delta= V2.^(1/2)-V1.^(1/2);  
        [~,q] = max(delta);     
        if mm~=q        
             aa(q)= aa(q) +1; %  YY(p,p)=Y(:,p)'*Y(:,p);
             aa(mm)= aa(mm) -1; %  YY(m,m)=Y(:,m)'*Y(:,m)
             ybby(mm)=V1(mm); % (30)
             ybby(q)=V2(q);
             Q(i,mm)=0;
             Q(i,q)=1;
             label(i)=q;
             SSQQ(:,mm)=SSQQ(:,mm)-SSS(:,i);% (29)
             SSQQ(:,q)=SSQQ(:,q)+SSS(:,i);     
        end  
     end   
%% compute objective function value
    obj(iter+1) = sum(ybby.^(1/2)) ;  
    if (iter>1 && abs(obj(iter)-obj(iter-1)) < 10^-6)
        break;
    end  
end 
    F=S*Q;
    [~,y]=max(F,[],2);
end

