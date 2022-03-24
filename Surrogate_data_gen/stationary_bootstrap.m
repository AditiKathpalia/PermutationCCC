function [x_surr]=stationary_bootstrap(x,numsurr)

%Generates stated number of surrogates of given time series using the
%stationary bootstrap method. For details of the method, check references: 
%1. Papana, A., Kyrtsou, C., Kugiumtzis, D., & Diks, C. (2017). 
%Assessment of resampling methods for causality testing: A note on the US 
%inflation behavior. PloS one, 12(7), e0180852.
%2. Politis, D. N., & Romano, J. P. (1994). The stationary bootstrap. 
%Journal of the American Statistical association, 89(428), 1303-1313.
%
% Aditi Kathpalia, NIAS


LEN=length(x);
p=0.1;

x_surr=zeros(numsurr,LEN);

for i=1:numsurr
    idx_1=randi([1,LEN],1,1);
    x_surr(i,1)=x(idx_1);
    idx=idx_1;
    for j=2:LEN
        rand_no=rand(1,1);
        if rand_no>p
            if idx<LEN
                idx=idx+1;
                x_surr(i,j)=x(idx);
            else
                idx=1;
                x_surr(i,j)=x(idx);
            end
        else
            idx=randi([1,LEN],1,1);
            x_surr(i,j)=x(idx);
        end
    end
    
end