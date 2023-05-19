clc;
clear all;
close all;

%% Loading Data and Model Identification

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

%%% Forming hk
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

% Form Hankel Tilde Matrix
for q=1:100
    p=1:2:N;
    index=1:100;
for r=1:100  
    H100_tilde(p(q):p(q)+1,p(r):p(r)+1)=h{r+index(q)};
end
end

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

%% Question 3: Calculating H_infinity
sys     = ss(A,B,C,D,ts);        % Define discrete state space system
H=tf(sys);                       % Transfer function of the state space
[ninf,fpeak] = hinfnorm(H);      % find inf norm using MATLAB build in function
[gam omeg]=infnorm_ds(A,B,C,D,10,0.1,0.0001) % find inf norm using function developed

%% Question 4
om = 2*pi*[0:(361)-1]/(ts*(361));  %%%% frequency vector in rad/s

% Emperical Freq Resp
y11f = fft(y11(mi:N))./fft(u1(mi:N));
y12f = fft(y12(mi:N))./fft(u2(mi:N));
y21f = fft(y21(mi:N))./fft(u1(mi:N));
y22f = fft(y22(mi:N))./fft(u2(mi:N));

% Model Freq Resp and SVD
for k=1:length(om)
temp1=C*inv(exp(i*ts*om(k))*eye(ns)-A)*B;
temp2=svd(temp1);
sv_freqres_m(k)=temp2(1,1);
temp3=svd([y11f(k) y12f(k); y21f(k) y22f(k)]);
sv_freqres_e(k)=temp3(1,1);
end

figure(1)
plot(om(1:182)/(2*pi),sv_freqres_m(1:182),'b*','LineWidth',1.2, 'MarkerSize',5)
hold on
plot(om(1:182)/(2*pi),sv_freqres_e(1:182),'r*','LineWidth',1.2, 'MarkerSize',5)
grid on
axis([0 20 0 0.5])
legend('Identified Model','Emperical','FontSize',12.5,'Interpreter','Latex')
xlabel('Frequency [Hz]','FontSize',12.5,'Interpreter','Latex')
ylabel('Singular Value','FontSize',12.5,'Interpreter','Latex')
title('Singular Values from 0 Hz to 20 Hz ','FontSize',12.5,'Interpreter','Latex')
axis square