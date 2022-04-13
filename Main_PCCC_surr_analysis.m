
close all;clear all;

% Load dataset. Demo rossler data is of dimension 2*32768*10
% First dimension is the variables, where, 1st row is independent and 2nd row is unidirectionally dependent on the first
% Second dimension is the length of time series (total =32768)
% Third dimension corresponds to the different realizations of the coupled
% system (total realizations=10)

load('Rossler_data_demo');

INFO_CCC=struct('L',25,'w',15,'delta',20);

% Compute time delay 'tau' using the first minimum of auto-mutual
% information or any other suitable method
% Compute embedding dimension
% 'emb_dim' using Cao's method or any other suitable method

emb_dim1=3;
tau1=5;            % Embedding dimension and delay of 1st time series

emb_dim2=3;
tau2=5;            % Embedding dimension and delay of 2nd time series

if tau1> tau2
    tau_to_use=tau1;
else
    tau_to_use=tau2;
end

n_sur=100;      % Number of surrogates for significance analysis
LEN=2048;        % Length of time series

if LEN-tau_to_use*emb_dim1> LEN-(INFO_CCC.w)
    nobs=LEN;
else
    nobs=LEN-tau_to_use*emb_dim1+(INFO_CCC.w);      % As only the 'cause' is to be embedded, remaining values in 'effect' can be used
end

INFO_CCC.N=nobs;                %CCC parameters

data(1,:)=ts_data(1,1:LEN,1);       % Loading first trial
data(2,:)=ts_data(2,1:LEN,1);

%Permutation binning the 'cause' time series, binning the 'effect' using equidistant with the same number of bins as in the permutation binned 'cause'
Y_CCC2=perm_binning_std(data(1,:),emb_dim1,tau1);   %Cause is embedded and permutation binned
X_CCC2=perm_binning_std(data(2,:),emb_dim2,tau2);

Y_CCC=Partition(data(1,:),factorial(emb_dim2));   %Effect is binned using equidistant binning
X_CCC=Partition(data(2,:),factorial(emb_dim1));

n_channels=2;
pointsToDiscard=0;

[CCC_value_dir1,k_right]=CCC_binned_seqs(X_CCC,Y_CCC2,INFO_CCC);        %To estimate CCC value from Y_CCC2 to X_CCC, call CCC function

[CCC_value_dir2,k_opp]=CCC_binned_seqs(Y_CCC,X_CCC2,INFO_CCC);              %To estimate CCC value from X_CCC2 to Y_CCC, call CCC function

% Surrogate data generation from both time series using AAFT/stationary
% bootstrap method, For rossler use AAFT, for real data bootstrap method
% [data1_s]=stationary_bootstrap(data(1,:),n_sur);
% [data2_s]=stationary_bootstrap(data(2,:),n_sur);

data1_s_temp = AAFTsur(data(1,:)',n_sur);
data2_s_temp = AAFTsur(data(2,:)',n_sur);
data1_s=data1_s_temp';
data2_s=data2_s_temp';
clear data1_s_temp data2_s_temp

parfor sur_no=1:n_sur
    
    Y_CCC2_s=perm_binning_std(data1_s(sur_no,:),emb_dim1,tau1);   % Permutation binning the cause surrogates
    X_CCC2_s=perm_binning_std(data2_s(sur_no,:),emb_dim2,tau2);
    
    Y_CCC_s=Partition(data1_s(sur_no,:),factorial(emb_dim2));   % Equidistant binning the effect surrogates
    X_CCC_s=Partition(data2_s(sur_no,:),factorial(emb_dim1));
    
    [CCC_value_dir1_sur(sur_no),k_right(sur_no)]=CCC_binned_seqs(X_CCC_s,Y_CCC2_s,INFO_CCC);        %To estimate CCC value, call CCC function
        
    [CCC_value_dir2_sur(sur_no),k_opp(sur_no)]=CCC_binned_seqs(Y_CCC_s,X_CCC2_s,INFO_CCC);
    
end

fprintf(1,'*****************Surrogate based significance testing result**************** \n\n');

%z-test based surrogate testing
mean_surr_CCC_dir1=mean(CCC_value_dir1_sur);
sigma_surr_CCC_dir1=std(CCC_value_dir1_sur);
[sig_dir1,p_val_dir1] = ztest(CCC_value_dir1,mean_surr_CCC_dir1,sigma_surr_CCC_dir1);

% Display significance testing result for direction 1
if sig_dir1==1
    fprintf(1,'Signicant PCCC in direction 1, i.e. from 1st variable/row of data to 2nd variable/row of data.\n\n');
else
    fprintf(1,'Insignicant PCCC in direction 1, i.e. from 1st variable/row of data to 2nd variable/row of data.\n\n');
end

fprintf(1,'****************** \n\n');

mean_surr_CCC_dir2=mean(CCC_value_dir2_sur);
sigma_surr_CCC_dir2=std(CCC_value_dir2_sur);
[sig_dir2,p_val_dir2] = ztest(CCC_value_dir2,mean_surr_CCC_dir2,sigma_surr_CCC_dir2);

% Display significance testing result for direction 2
if sig_dir2==1
    fprintf(1,'Signicant PCCC in direction 2, i.e. from 2nd variable/row of data to 1st variable/row of data.\n\n');
else
    fprintf(1,'Insignicant PCCC in direction 2, i.e. from 2nd variable/row of data to 1st variable/row of data.\n\n');
end

fprintf(1,'***************************************************************************** \n');
