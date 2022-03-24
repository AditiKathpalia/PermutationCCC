function [New_symb_seq]=perm_binning_std(ts,n,tau)
% Transforms given one dimensional time series (coming from an
% n-dimensional dynamical system) to a one-dimensional ordinal pattern
% sequence
% Inputs -
% ts - Input real-valued time series
% n - embedding dimension
% tau - embedding delay
%
% Outputs -
% New_symb_seq - One dimensional symbolic sequence with ordinal patterns (represented by n! symbols)
%
% Aditi Kathpalia, ICS, CAS, May 2021, Last modified: March 2022

N=length(ts);
symb_x=zeros(1,n);

for i=1:N-n*tau+1
    x=ts((i:tau:i+n*tau-1));    % Taking n values from ts at a step-size tau to form vector 'x', Moving window step for vectors =1
    [A,I]=sort(x);               % Sort values in x
    for k=1:n
        symb_x(I(k))=k-1;          % Label lowest value with lowest integer 'k'
    end
    
    ts_m(i,:)=symb_x;               % Ordinal patterns at time step 'i'
    d(i) = bi2de(ts_m(i,:),n,'left-msb');   % Converts the base-n vector to a decimal integer (1D conversion of vector)
    
end

d_unique=unique(d);         % Find the unique integers/ ordinal patterns (in ascending order)

Symbols=1:factorial(n);     % Define symbols to represent the final sequence

for k=1:length(d)
    idx=find(d_unique==d(k));       % find the index of d(k) in d_unique
    d_new(k)=Symbols(idx);          % Replace with desired symbol value at the same index
end

New_symb_seq=d_new;

