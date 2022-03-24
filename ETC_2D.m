function [N, EntLen_array, INFO1, INFO2] = ETC_2D(InputTS1,InputTS2,NumBins)
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
%  >> [N, EntLen_array] = ETC([1 2 1 2 2 2 2 ],0)
%  gives :  N=5 ,  EntLen_array = [6.0418    4.8548    6.0000    2.7549    2.0000  0]
%
% Example2:
% >> [N, EntLen_array] = ETC([0:0.1:0.99],2)
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
% v3.0 (Feb 10, 2016) - Nithin@NIAS, 2D version

N = -1; %ETC measure using NSRPS is initialized to -1
L = length(InputTS1);  % should be same as length(InputTS2)


if NumBins~=0
    [SymSeq1, INFO1] = Partition(InputTS1,NumBins); % input is a time series of real values
    [SymSeq2, INFO2] = Partition(InputTS2,NumBins); % input is a time series of real values
else
    SymSeq1=InputTS1;  % input is symbolic sequence
    INFO1.NumSymbols = length(unique(SymSeq1));
    SymSeq2=InputTS2;  % input is symbolic sequence
    INFO2.NumSymbols = length(unique(SymSeq2));
end
% % % to get rid of zeros in the symbolic sequence, if any
% % minY = min(SymSeq); 
% % y = SymSeq - minY;
% % y = y+1; 
% % SymSeq=y;

if ~isempty(find(SymSeq1==0, 1))
    fprintf(1,'\n ERROR, 0 found \n');
end

if ~isempty(find(SymSeq2==0, 1))
    fprintf(1, '\n ERROR, 0 found \n');
end

% The main loop for NSRPS iteration
N=0; % ETC measure 'N'
Hnew = ShannonEntropy2(SymSeq1, SymSeq2); % Shannon entropy of the symbolic sequence 
Len = length(SymSeq1);
EntLen_array=[Hnew*Len];

while (Hnew > 1e-6) && (Len>1) 
[Pair1, Pair2] = FindPair2(SymSeq1, SymSeq2);  % find the pair of symbols with maximum frequency

[SymSeqNew1, SymSeqNew2] = Substitute2(SymSeq1, SymSeq2, Pair1, Pair2);

% [SymSeqNew1] = Substitute(SymSeq1, Pair1);  % substitute the pair with a new symbol
% [SymSeqNew2] = Substitute(SymSeq2, Pair2);  % substitute the pair with a new symbol


Hnew = ShannonEntropy2(SymSeqNew1, SymSeqNew2); % Shannon entropy of the new sequence 
Len=length(SymSeqNew1); 
EntLen_array = [EntLen_array Hnew*Len];
N=N+1; % ETC measure incremented by 1
SymSeq1 = SymSeqNew1;
SymSeq2 = SymSeqNew2;
clear SymSeqNew1 SymSeqNew1
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
function [Pair1, Pair2] = FindPair2(SymSeq1, SymSeq2)
% Computes the pair with maximum frequency in the symbolic seq.
L1 = length(SymSeq1); % assume that both symbolic sequence are of same length
L2 = length(SymSeq2); % assume that both symbolic sequence are of same length

if L1 ~= L2
    disp('Error, mismatch in lengths');
    return;
end

L = L1;
X = SymSeq1;
Y = SymSeq2;
s1 = struct('pair2d', [X(1) X(2) Y(1) Y(2)], 'count', 1);
S{1} = s1;
INDX = 1;
for i = 2:L-1
    pair_buf = [X(i) X(i+1) Y(i) Y(i+1)];
    
    for j = 1:length(S)
        RES = isequal(S{j}.pair2d,pair_buf);
        if RES==1;
           S{j}.count = S{j}.count+1; 
        end
    end
    
    if RES==0
        s_buf = struct('pair2d', pair_buf, 'count', 1);
        S{INDX+1} = s_buf;
        INDX = INDX+1;
        clear s_buf
    end
    
    clear pair_buf RES
end

% Now find maximum count 4-tuple
s_buf = S{1};
MAX_count = s_buf.count;
MAX_pair = s_buf.pair2d;
for k = 2:length(S)

    if S{k}.count > MAX_count
        s_buf = S{k};
        MAX_count = s_buf.count;
        MAX_pair = s_buf.pair2d;
    end
end

Pair1 = [MAX_pair(1) MAX_pair(2)];
Pair2 = [MAX_pair(3) MAX_pair(4)];

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
function [SymSeqNew1, SymSeqNew2] = Substitute2(SymSeq1, SymSeq2, Pair1, Pair2) 
% Pair substitution step of NSRPS - 2D version

SymSeqNew1 = [];
SymSeqNew2 = [];

L = length(SymSeq1);
Alphabet = unique([SymSeq1, SymSeq2]);
M = length(Alphabet);
I = max(Alphabet);
RepSym = I+1; % New Symbol
indx=1;

while (indx<L)
    a1=SymSeq1(indx);  a2=SymSeq2(indx);
    b1=SymSeq1(indx+1); b2=SymSeq2(indx+1);
    
    
    if (a1==Pair1(1))&&(b1==Pair1(2)) && (a2==Pair2(1))&&(b2==Pair2(2))
        SymSeqNew1 = [SymSeqNew1 RepSym];
        SymSeqNew2 = [SymSeqNew2 RepSym];
        indx=indx+1;
    else
        SymSeqNew1=[SymSeqNew1 a1];
        SymSeqNew2=[SymSeqNew2 a2];
    end
    indx=indx+1;
    if indx==L;
        SymSeqNew1 =[SymSeqNew1 SymSeq1(end)];
        SymSeqNew2 =[SymSeqNew2 SymSeq2(end)];
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

%----------------------------------------------------
function [H] = ShannonEntropy2(SymSeq1, SymSeq2)
% Computes the Shannon Entropy of a pair of symbolic sequences (values need
% to be strictly > 0, else it will give an error)
% also, the two symbolic sequences need to have the same length


% % to make enteries in the symbolic sequence stictly >0
% minY = min(SymSeq); 
% y = SymSeq - minY;
% y = y+1; 

y1 = SymSeq1; % has to be integer entries > 0
y2 = SymSeq2; % has to be integer entries > 0
H = 0;
L = length(y1);
Num1 = max(y1); Num2 = max(y2);
prob2D = zeros(Num1, Num2);

for i = 1:L;
    prob2D(y1(i),y2(i))= prob2D(y1(i),y2(i))+1;
end

prob2D = prob2D/L;

for i1 = 1:Num1;
    for i2 = 1:Num2;
        if prob2D(i1,i2)~=0;
            H = H - prob2D(i1,i2)*log2(prob2D(i1,i2));
        end
    end
end
end
%----------------------------------------------------


% %----------------------------------------------------
% function [Pair1, Pair2] = FindPair2(SymSeq1, SymSeq2, NumBins)
% % Computes the pair with maximum frequency in the symbolic seq.
% Alphabet1 = unique(SymSeq1);
% Alphabet2 = unique(SymSeq2);
% M1 = max(Alphabet1);
% M2 = max(Alphabet2);
% M = NumBins; %max(M1,M2);
% Count_Array = zeros(1,M^4);
% 
% A = SymSeq1-1;
% B = SymSeq2-1;
% 
% L = length(SymSeq1);
% indx = 1;
% while (indx < L);
%     a = [A(indx) A(indx+1)];
%     b = [B(indx) B(indx+1)];
%     
%     VAL = a(1) + a(2)*NumBins + b(1)*NumBins^2 + b(2)*NumBins^3;
%     Count_Array(VAL+1)  = Count_Array(VAL+1) + 1; % to take are that VAL can be zero
%     indx = indx+1;
% end
% 
% [m,INDEX]=max(Count_Array(:));
% 
% PAIR = [];
% X = INDEX-1;
% for i = 1:4
%    PAIR = [PAIR mod(X,NumBins)];
%    X = floor(X/NumBins);
% end
% 
% Pair1 = [PAIR(1) PAIR(2)];
% Pair2 = [PAIR(3) PAIR(4)];
% 
% Pair1 = Pair1+1;
% Pair2 = Pair2+1;
% 
% end
% %----------------------------------------------------

