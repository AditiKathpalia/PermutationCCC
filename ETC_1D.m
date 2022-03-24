function [N, EntLen_array, INFO] = ETC_1D(InputTS,NumBins)
% Computation of "Effort-To-Compress" (ETC) measure using Non-Sequential Recursive Pair
% Substitution (NSRPS)
%
%   INPUTS: 
%   InputTS - input time series of real numbers or a symbolic sequence of
%             integers
%   NumBins - Number of Bins to partition the input time series to create a symbolic sequence,
%             if this is set to 0 then input is already a symbolic sequence
%
%   OUTPUTS:
%   N - ETC measure = number of iterations for entropy -> 0 in the NSRPS
%                     algorithm
%   EntLen_array - Array of entropy x length values for the NSRPS algorithm    
%
%  Example1: 
%  >> [N, EntLen_array] = ETC_1D([1 2 1 2 2 2 2 ],0)
%  gives :  N=5 ,  EntLen_array = [6.0418    4.8548    6.0000    2.7549    2.0000  0]
%
% Example2:
% >> [N, EntLen_array] = ETC_1D([0:0.1:0.99],2)
% gives: N=7 , EntLen_array =[10.0000   10.3904   11.5098    9.6096    8.0000    4.7549    2.0000 0]
% 
% 
% by: Nithin Nagaraj, Feb. 2013
% Amrita Vishwa Vidyapeetham, Amritapuri, Kollam, Kerala, INDIA
%
% Ref: 1. Nithin Nagaraj, Karthi Balasubramanian, Sutirth Dey, A New Complexity Measure For Time Series Analysis and
% Classification, Eur. Phys. J. Special Topics 222, 847–860 (2013)
%
% 2. W. Ebeling, M.A. Jim\'{e}nez-Monta\~{n}o, On Grammars, Complexity, and
% Information Measures of Biological Macromolecules, Mathematical Biosciences 52, (1980) 53-71.
%
% The latest version of this program can be downloaded from the following websites:
%
% URL1: https://sites.google.com/site/nithinnagaraj2/journal/etc
% URL2: https://sites.google.com/a/acads.iiserpune.ac.in/sdlab/publications?pli=1
%
% This program may be used for research purposes only. We provide no
% guarantees on the validity of the output, or the interpretation of
% results. Please do not remove this notice. 
%
% v1.0 (Feb 15, 2013)
% v1.1 (May 31, 2013) 
% v2.0 (Dec 29, 2015) - Nithin@NIAS, added INFO structure

N = -1; %ETC measure using NSRPS is initialized to -1
L = length(InputTS);


if NumBins~=0
    [SymSeq, INFO] = Partition(InputTS,NumBins); % input is a time series of real values
else
    SymSeq=InputTS;  % input is symbolic sequence
    INFO.NumSymbols = length(unique(SymSeq));
end

% % to get rid of zeros in the symbolic sequence, if any
% minY = min(SymSeq); 
% y = SymSeq - minY;
% y = y+1; 
% SymSeq=y;

if ~isempty(find(SymSeq==0, 1))
   fprintf(1,'\n ERROR, 0 found \n');
end

% The main loop for NSRPS iteration
N=0; % ETC measure 'N'
Hnew = ShannonEntropy(SymSeq); % Shannon entropy of the symbolic sequence 
Len = length(SymSeq);
EntLen_array=[Hnew*Len];

while (Hnew > 1e-6) && (Len>1) 
Pair = FindPair(SymSeq);  % find the pair of symbols with maximum frequency
SymSeqNew = Substitute(SymSeq,Pair);  % substitute the pair with a new symbol
Hnew = ShannonEntropy(SymSeqNew); % Shannon entropy of the new sequence 
Len=length(SymSeqNew); 
EntLen_array = [EntLen_array Hnew*Len];
N=N+1; % ETC measure incremented by 1
SymSeq = SymSeqNew;
clear SymSeqNew;
end

end
%--------------end of ETC function-------------------

%----------------------------------------------------
% Subroutine functions called by ETC function
%----------------------------------------------------
function [SymSeq, INFO] = Partition(InputTS,NumBins)
% Coverts an input time series into a symbolic sequence 

x=InputTS; SymSeq=x;
L = length(InputTS);
Range = max(x(:))+1e-6-min(x(:));
Delta = Range/NumBins;
x = x - min(x);
for i = 1:L;
    SymSeq(i) = floor(x(i)/Delta);
end
SymSeq = SymSeq+1;

%>>>>>>>>> Dec 29 Nithin@NIAS INFO structure added
%This is for bookeeping - Information about max, min of data-set or InputTS
INFO.InputTS = InputTS;
INFO.NumBins = NumBins;
INFO.min = min(InputTS(:));
INFO.max = max(InputTS(:));  
INFO.range = Range;
INFO.delta = Delta;
INFO.Partitions = {}; % goes like [min Delta], (Delta, 2Delta], (2Delta, 3Delta]....((NumBins-1)*Delta ,Max]
INFO.SymSeq = SymSeq;

INFO.Partitions = ['[', num2str(INFO.min), ',', num2str(INFO.delta), '],'];
for i = 2:NumBins-1;
    INFO.Partitions = [INFO.Partitions, ['(', num2str((INFO.delta)*(i-1)), ',', num2str((INFO.delta)*(i)), '],'] ];
end
 INFO.Partitions = [INFO.Partitions, ['(', num2str((INFO.delta)*(NumBins-1)), ',', num2str(INFO.max), ']'] ];
%>>>>>>>>>

end
%----------------------------------------------------


%----------------------------------------------------
function [Pair] = FindPair(SymSeq)
% Computes the pair with maximum frequency in the symbolic seq.

Pair = [];
Alphabet = unique(SymSeq);
M = max(Alphabet);
Count_Array = zeros(M,M);
L = length(SymSeq);
indx=1;
while (indx<L)
    a=SymSeq(indx);
    b=SymSeq(indx+1);
    Count_Array(a,b)=Count_Array(a,b)+1;
    if a==b
        if indx<L-1
            if (SymSeq(indx+2)==a)
                indx=indx+1;
            end
        end
    end
    indx=indx+1;
end
[m,indx]=max(Count_Array(:));
Pair1=mod((indx-1),M)+1;
Pair2=floor((indx-1)/M)+1;
Pair = [Pair1 Pair2];
end
%----------------------------------------------------

%----------------------------------------------------
function [SymSeqNew]= Substitute(SymSeq,Pair)
% Pair substitution step of NSRPS

SymSeqNew = [];
L = length(SymSeq);
Alphabet = unique(SymSeq);
M = length(Alphabet);
I = max(Alphabet);
RepSym = I+1; % New Symbol
indx=1;
while (indx<L)
    a=SymSeq(indx);
    b=SymSeq(indx+1);
    if (a==Pair(1))&&(b==Pair(2))
        SymSeqNew = [SymSeqNew RepSym];
        indx=indx+1;
    else
        SymSeqNew=[SymSeqNew a];
    end
    indx=indx+1;
    if indx==L;
        SymSeqNew =[SymSeqNew SymSeq(end)];
    end 
end
end
%----------------------------------------------------

%----------------------------------------------------
function [H] = ShannonEntropy(SymSeq)
% Computes the Shannon Entropy of a symbolic sequence

% to make enteries in the symbolic sequence stictly >0
minY = min(SymSeq); 
y = SymSeq - minY;
y = y+1; 

y=SymSeq;
H = 0;
L = length(y);
Num = max(y); 
prob = zeros(1,Num);

for i = 1:L;
    prob(y(i))= prob(y(i))+1;
end

prob = prob/L;

for i = 1:Num;
    if prob(i)~=0;
        H = H - prob(i)*log2(prob(i));
    end
end
end
%----------------------------------------------------