function[CCC_value,k]=CCC_binned_seqs(X,Y,INFO_CCC)
%Estimates causality using Compression Complexity Causality(CCC) measure
%from input time series Y to X when the time series are already binned
%Outputs-
%   1. CCC_value - Estimated causality using Compression Complexity Causality(CCC)
%   2. k - no. of times the difference between complexity rate based on its own past and joint is estimated from the given series
%
%Inputs -
%   1. Binned sequences X, Y to check causation from Y to X
%   2. INFO_CCC is a structure taking in input parameters -
%      N=length of input sequences, L=Past window size, w=current window size, delta=Jump step size for each window
%Calls functions ETC_1D, ETC_2D. Check the readme file for more details.
%
%
% Aditi Kathpalia, ICS, CAS, May 2021, Last modified: March 2022

LEN_past=INFO_CCC.L;
ADD_meas=INFO_CCC.w;
Num_bins=0;             %Set the number of bins required for ETC as zero as the sequences are already binned
LEN=INFO_CCC.N;
STEP_size=INFO_CCC.delta;

LEN_to_check=LEN_past+ADD_meas;
sequence1=X;
sequence2=Y;

% Compute dynamical complexities and their differences

k=1;
for i=1:STEP_size:LEN-LEN_to_check
    
    ETC_ini(k)=ETC_1D(sequence1(i:i+LEN_past-1),Num_bins);
    ETC_ini(k)=ETC_ini(k)/(LEN_past-1);
    ETC_ini_2D(k)=ETC_2D(sequence1(i:i+LEN_past-1), sequence2(i:i+LEN_past-1),Num_bins);
    ETC_ini_2D(k)=ETC_ini_2D(k)/(LEN_past-1);
    ETC_fin(k)=ETC_1D(sequence1(i:i+LEN_to_check-1),Num_bins);
    ETC_fin(k)=ETC_fin(k)/(LEN_to_check-1);
    ETC_fin_2D(k)=ETC_2D(sequence1(i:i+LEN_to_check-1), [sequence2(i:i+LEN_past-1) sequence1(i+LEN_past:i+LEN_to_check-1)],Num_bins);
    ETC_fin_2D(k)=ETC_fin_2D(k)/(LEN_to_check-1);
    
    delta_ETC(k)=ETC_fin(k)-ETC_ini(k);
    delta_ETC_2D(k)=ETC_fin_2D(k)-ETC_ini_2D(k);
    
    k=k+1;
end
%      figure;plot(delta_ETC,'b');
%     hold on;plot(delta_ETC_2D,'r');

% To estimate average difference between dynamical complexities and hence CCC

sum_delta_ETC=sum(delta_ETC);
sum_delta_ETC_2D=sum(delta_ETC_2D);

CCC_value=(sum_delta_ETC-sum_delta_ETC_2D)/(k-1);
