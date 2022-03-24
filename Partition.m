function [ y ] = Partition( x, Num_bins)
% Time series x is transformed to a symbolic sequence y using no. of bins ='Num_bins'


if min(x)==max(x)
    y=x;
else
    x1=x-min(x);
    x1=x1./(max(x)-min(x));
    
    y= ceil(x1.*Num_bins);
end

y(find(y==0))=1;


end

