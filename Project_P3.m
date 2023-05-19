clear all;
close all;
clc

%% Loading Data from White Noise Response
load u_rand.mat
y1 = 0.5*u_rand.Y(3).Data;
y2 = 0.5*u_rand.Y(4).Data;
u1 = 0.5*u_rand.Y(1).Data;
u2 = 0.5*u_rand.Y(2).Data;
ts = 1/40;
N = length(y1);
t = [0:N-1]*ts - 1;


%% Ques #1

% Define var for y auto correlation
yy11=zeros(1,1);
yy21=zeros(1,1);
yy12=zeros(1,1);
yy22=zeros(1,1);

p=12000; % Half of data point of y
  
% Calculating for y_RMS using y data points
    for k=-p:p
        
        temp=[y1(k+12001) ; y2(k+12001)]*[y1(k+12001) y2(k+12001)];
        yy11=yy11+temp(1,1); 
        yy21=yy21+temp(2,1);
        yy12=yy12+temp(1,2);
        yy22=yy22+temp(2,2);

    end

% y_RMS 
y_RMS=sqrt(trace((1/(2*p))*[yy11 yy12; yy21 yy22])) 

% Calculate y_RMS using R_yy[0]
% Define var for R_yy
Ryy11=zeros(1,401);
Ryy21=zeros(1,401);
Ryy12=zeros(1,401);
Ryy22=zeros(1,401);

for k=-200:200 % calculate from t=-5 to t=5
    
    for q=-11800:11800
        
        temp=[y1(k+q+12001) ; y2(k+q+12001)]*[y1(q+12001) y2(q+12001)];
        Ryy11(k+201)=Ryy11(k+201)+temp(1,1); 
        Ryy21(k+201)=Ryy21(k+201)+temp(2,1);
        Ryy12(k+201)=Ryy12(k+201)+temp(1,2);
        Ryy22(k+201)=Ryy22(k+201)+temp(2,2);

    end
    Ryy11(k+201)=(1/(2*11800+1))*Ryy11(k+201); 
    Ryy21(k+201)=(1/(2*11800+1))*Ryy21(k+201);
    Ryy12(k+201)=(1/(2*11800+1))*Ryy12(k+201);
    Ryy22(k+201)=(1/(2*11800+1))*Ryy22(k+201);
end

% y_RMS from R_yy[0]
y_RMS_Ryy=sqrt(trace([Ryy11(201) Ryy12(201);Ryy21(201) Ryy22(201)]))


%% Load Data from Impulse Response / Ques #2

load u1_impulse.mat

y11 = u1_impulse.Y(3).Data;
y21 = u1_impulse.Y(4).Data;
u1 = u1_impulse.Y(1).Data;  %%% note that the pulse magnitude is 5

[m,mi] = max(u1>0);  %%% find index where pulse occurs

load u2_impulse.mat

y12 = u2_impulse.Y(3).Data;
y22 = u2_impulse.Y(4).Data;
u2 = u2_impulse.Y(2).Data;

%%% remove any offsets in output data using data prior to pulse application
y11 = y11 - mean(y11([1:mi-1]));
y12 = y12 - mean(y12([1:mi-1]));
y21 = y21 - mean(y21([1:mi-1]));
y22 = y22 - mean(y22([1:mi-1]));

%%% rescale IO data so that impulse input has magnitude 1
y11 = y11/max(u1);
y12 = y12/max(u2);
y21 = y21/max(u1);
y22 = y22/max(u2);
u1 = u1/max(u1);
u2 = u2/max(u2);

ts = 1/40;  %%%% sample period
N = length(y11);  %%%% length of data sets
t = [0:N-1]*ts - 1;

% Forming hk
for k=1:N-mi
h{k}=[y11(k+mi) y12(k+mi);y21(k+mi) y22(k+mi)];
end

% Form Hankel Matrix
for q=1:100
    p=1:2:N;
    index=0:99;
for r=1:100  
    H100(p(q):p(q)+1,p(r):p(r)+1)=h{r+index(q)};
end
end

% Form Hankel tilde Matrix
for q=1:100
    p=1:2:N;
    index=1:100;
for r=1:100  
    H100_tilde(p(q):p(q)+1,p(r):p(r)+1)=h{r+index(q)};
end
end

% Find A,B & C

m=2;                    % Number of inputs
q=2;                    % Number of Output    
ns=7;                   % System State
n=100;                  % Size of Hankel Matrix

[U,sig,V]=svd(H100);    % SV Decomposition of Hankel
  
On=U(1:m*n,1:ns)*sig(1:ns,1:ns);     % Observability Mat
Cn=(V(1:n*q,1:ns))';                 % Controllability Mat

C=On(1:m,:);                         % C from Observability Mat
B=Cn(:,1:q);                         % B from Controllability Mat

li_On=inv(On'*On)*On';              % Left-inverse of Observability Mat
ri_Cn=Cn'*inv(Cn*Cn');              % Right-inverse of Controllability Mat

A=li_On*H100_tilde*ri_Cn;           % Calculating A
D=zeros(2,2);                       % Define D

%% Ques #2
sys     = ss(A,B,C,D,ts);   % Define discrete state space system
Gc      = gram(sys,'c');    % observability mat
Go      = gram(sys,'o');    % controllability mat
PH2_a   = sqrt(trace(B'*Go*B))  % Using eqn 9 with observability mat
PH2_b   = sqrt(trace(C*Gc*C'))  % Using eqn 10 with observability mat


%% Question 3: H2 using pulse response

PH2_impulse=0; % Define variable for the euclidean norm of the system

for k=1:length(h)
    
    % Calclate frobenious norm for each time step of impulse resp matrix
    temp=norm(h{k},'fro')^2; 
    PH2_impulse=PH2_impulse+temp; 
    
end

PH2_impulse=sqrt(PH2_impulse)
