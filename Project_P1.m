clear all;
close all;
clc

%% Loading Data
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

%% Plotting Impulse Response
figure(1);
subplot(311)
plot(t,u1,'b*','LineWidth',2)
ylabel('$u_1$ (volts)','FontSize',14,'Interpreter','Latex');
grid on
axis([-0.2 2 -0.1 1.1])

subplot(312)
plot(t,y11,'r*','LineWidth',1)
ylabel('$y_1$ (volts)','FontSize',14,'Interpreter','Latex');
grid on
axis([-0.2 2 -0.1 0.1])

subplot(313)
plot(t,y21,'r*','LineWidth',1)
ylabel('$y_2$ (volts)','FontSize',14,'Interpreter','Latex');
grid on
xlabel('second','FontSize',14)
set(gca,'FontSize',14)
axis([-0.2 2 -0.1 0.1])

figure(2);
subplot(311)
plot(t,u2,'b*','LineWidth',2)
ylabel('$u_2$ (volts)','FontSize',14,'Interpreter','Latex');
grid on
axis([-0.2 2 -0.1 1.1])

subplot(312)
plot(t,y12,'r*','LineWidth',1)
ylabel('$y_1$ (volts)','FontSize',14,'Interpreter','Latex');
grid on
axis([-0.2 2 -0.1 0.1])

subplot(313)
plot(t,y22,'r*','LineWidth',1)
ylabel('$y_2$ (volts)','FontSize',14,'Interpreter','Latex');
grid on
xlabel('second','FontSize',14)
set(gca,'FontSize',14)
axis([-0.2 2 -0.1 0.1])

%% Form Hankel Matrix / Ques #1

for k=1:N-mi
h{k}=[y11(k+mi) y12(k+mi);y21(k+mi) y22(k+mi)]; % Forming hk or impulse mat
end

%%% Hankel Matrix size=100
for q=1:100
    p=1:2:N;
    index=0:99;
for r=1:100  
    H100(p(q):p(q)+1,p(r):p(r)+1)=h{r+index(q)};
end
end

%%% Hankel Tilde Matrix size=100
for q=1:100
    p=1:2:N;
    index=1:100;
for r=1:100  
    H100_tilde(p(q):p(q)+1,p(r):p(r)+1)=h{r+index(q)};
end
end

%% Plot SV of H100 / Ques #1
svdy=svd(H100);              %Finding the SV
svdx=linspace(1,200,200)';   

figure(3)
semilogy(svdx,svdy,'o', 'MarkerFaceColor', 'b')
grid on
axis([0 40 1e-3 1])
xlabel('Singular Value Index')
ylabel('Hankel Singular Value')
title('H_{100} Singular Values for the First 40')

%% Find A,B & C / Ques #1
m=2;                    % Number of inputs
q=2;                    % Number of Output    
ns=[6 7 10 20 8];       % System State
n=100;                  % Size of Hankel Matrix

[U,sig,V]=svd(H100);    % SV Decomposition of Hankel

for y=1:length(ns)
    
On{y}=U(1:m*n,1:ns(y))*sig(1:ns(y),1:ns(y));    % Observability Mat
Cn{y}=(V(1:n*q,1:ns(y)))';                      % Controllability Mat

C{y}=On{y}(1:m,:);  % C from Observability Mat
B{y}=Cn{y}(:,1:q);  % B from Controllability Mat

li_On{y}=inv(On{y}'*On{y})*On{y}'; % Left-inverse of Observability Mat
ri_Cn{y}=Cn{y}'*inv(Cn{y}*Cn{y}'); % Right-inverse of Controllability Mat

A{y}=li_On{y}*H100_tilde*ri_Cn{y}; % Calculating A
end

max(abs(eig(A{3}))) % Check for stability (must less than 1)

%% Impulse Response from State Space / Ques #2

for k=1:N-mi
impul{k}=C{1}*(A{1}^(k-1))*B{1};    % Impulse for ns=6
h11a(k) = impul{k}(1,1);
h12a(k) = impul{k}(1,2);
h21a(k) = impul{k}(2,1);
h22a(k) = impul{k}(2,2);
end

for k=1:N-mi
impul{k}=C{2}*(A{2}^(k-1))*B{2};     % Impulse for ns=7
h11b(k) = impul{k}(1,1);
h12b(k) = impul{k}(1,2);
h21b(k) = impul{k}(2,1);
h22b(k) = impul{k}(2,2);
end

for k=1:N-mi
impul{k}=C{3}*(A{3}^(k-1))*B{3};     % Impulse for ns=10
h11c(k) = impul{k}(1,1);
h12c(k) = impul{k}(1,2);
h21c(k) = impul{k}(2,1);
h22c(k) = impul{k}(2,2);
end

for k=1:N-mi
impul{k}=C{4}*(A{4}^(k-1))*B{4};     % Impulse for ns=20
h11d(k) = impul{k}(1,1);
h12d(k) = impul{k}(1,2);
h21d(k) = impul{k}(2,1);
h22d(k) = impul{k}(2,2);
end

% Impulse for ns=6
figure(4)
subplot(411)
plot(t(42:N),h11a,'b*','LineWidth',1)
hold on
plot(t(42:N),y11(42:N),'r*','LineWidth',1)
ylabel('$y_{11}$ (volts)','FontSize',14,'Interpreter','Latex');
legend('Calculated h(t)','Measured h(t)')
grid on
xlabel('second','FontSize',14)
set(gca,'FontSize',14)
axis([0 2 -0.1 0.1])
hold off


subplot(412)
plot(t(42:N),h12a,'b*','LineWidth',1)
hold on
plot(t(42:N),y12(42:N),'r*','LineWidth',1)
ylabel('$y_{12}$ (volts)','FontSize',14,'Interpreter','Latex');
legend('Calculated h(t)','Measured h(t)')
grid on
xlabel('second','FontSize',14)
set(gca,'FontSize',14)
axis([0 2 -0.1 0.1])
hold off

subplot(413)
plot(t(42:N),h21a,'b*','LineWidth',1)
hold on
plot(t(42:N),y21(42:N),'r*','LineWidth',1)
ylabel('$y_{21}$ (volts)','FontSize',14,'Interpreter','Latex');
legend('Calculated h(t)','Measured h(t)')
grid on
xlabel('second','FontSize',14)
set(gca,'FontSize',14)
axis([0 2 -0.1 0.1])
hold off

subplot(414)
plot(t(42:N),h22a,'b*','LineWidth',1)
hold on
plot(t(42:N),y22(42:N),'r*','LineWidth',1)
ylabel('$y_{22}$ (volts)','FontSize',14,'Interpreter','Latex');
legend('Calculated h(t)','Measured h(t)')
grid on
xlabel('second','FontSize',14)
set(gca,'FontSize',14)
axis([0 2 -0.1 0.1])
hold off

% Impulse for ns=7
figure(5)
subplot(411)
plot(t(42:N),h11b,'b*','LineWidth',1)
hold on
plot(t(42:N),y11(42:N),'r*','LineWidth',1)
ylabel('$y_{11}$ (volts)','FontSize',14,'Interpreter','Latex');
legend('Calculated h(t)','Measured h(t)')
grid on
xlabel('second','FontSize',14)
set(gca,'FontSize',14)
axis([0 2 -0.1 0.1])
hold off

subplot(412)
plot(t(42:N),h12b,'b*','LineWidth',1)
hold on
plot(t(42:N),y12(42:N),'r*','LineWidth',1)
ylabel('$y_{12}$ (volts)','FontSize',14,'Interpreter','Latex');
legend('Calculated h(t)','Measured h(t)')
grid on
xlabel('second','FontSize',14)
set(gca,'FontSize',14)
axis([0 2 -0.1 0.1])
hold off

subplot(413)
plot(t(42:N),h21b,'b*','LineWidth',1)
hold on
plot(t(42:N),y21(42:N),'r*','LineWidth',1)
ylabel('$y_{21}$ (volts)','FontSize',14,'Interpreter','Latex');
legend('Calculated h(t)','Measured h(t)')
grid on
xlabel('second','FontSize',14)
set(gca,'FontSize',14)
axis([0 2 -0.1 0.1])
hold off

subplot(414)
plot(t(42:N),h22b,'b*','LineWidth',1)
hold on
plot(t(42:N),y22(42:N),'r*','LineWidth',1)
ylabel('$y_{22}$ (volts)','FontSize',14,'Interpreter','Latex');
legend('Calculated h(t)','Measured h(t)')
grid on
xlabel('second','FontSize',14)
set(gca,'FontSize',14)
axis([0 2 -0.1 0.1])
hold off

% Impulse for ns=10
figure(6)
subplot(411)
plot(t(42:N),h11c,'b*','LineWidth',1)
hold on
plot(t(42:N),y11(42:N),'r*','LineWidth',1)
ylabel('$y_{11}$ (volts)','FontSize',14,'Interpreter','Latex');
legend('Calculated h(t)','Measured h(t)')
grid on
xlabel('second','FontSize',14)
set(gca,'FontSize',14)
axis([0 2 -0.1 0.1])
hold off

subplot(412)
plot(t(42:N),h12c,'b*','LineWidth',1)
hold on
plot(t(42:N),y12(42:N),'r*','LineWidth',1)
ylabel('$y_{12}$ (volts)','FontSize',14,'Interpreter','Latex');
legend('Calculated h(t)','Measured h(t)')
grid on
xlabel('second','FontSize',14)
set(gca,'FontSize',14)
axis([0 2 -0.1 0.1])
hold off

subplot(413)
plot(t(42:N),h21c,'b*','LineWidth',1)
hold on
plot(t(42:N),y21(42:N),'r*','LineWidth',1)
ylabel('$y_{21}$ (volts)','FontSize',14,'Interpreter','Latex');
legend('Calculated h(t)','Measured h(t)')
grid on
xlabel('second','FontSize',14)
set(gca,'FontSize',14)
axis([0 2 -0.1 0.1])
hold off

subplot(414)
plot(t(42:N),h22c,'b*','LineWidth',1)
hold on
plot(t(42:N),y22(42:N),'r*','LineWidth',1)
ylabel('$y_{22}$ (volts)','FontSize',14,'Interpreter','Latex');
legend('Calculated h(t)','Measured h(t)')
grid on
xlabel('second','FontSize',14)
set(gca,'FontSize',14)
axis([0 2 -0.1 0.1])
hold off

% Impulse for ns=20
figure(7)
subplot(411)
plot(t(42:N),h11d,'b*','LineWidth',1)
hold on
plot(t(42:N),y11(42:N),'r*','LineWidth',1)
ylabel('$y_{11}$ (volts)','FontSize',14,'Interpreter','Latex');
legend('Calculated h(t)','Measured h(t)')
grid on
xlabel('second','FontSize',14)
set(gca,'FontSize',14)
axis([0 2 -0.1 0.1])
hold off

subplot(412)
plot(t(42:N),h12d,'b*','LineWidth',1)
hold on
plot(t(42:N),y12(42:N),'r*','LineWidth',1)
ylabel('$y_{12}$ (volts)','FontSize',14,'Interpreter','Latex');
legend('Calculated h(t)','Measured h(t)')
grid on
xlabel('second','FontSize',14)
set(gca,'FontSize',14)
axis([0 2 -0.1 0.1])
hold off

subplot(413)
plot(t(42:N),h21d,'b*','LineWidth',1)
hold on
plot(t(42:N),y21(42:N),'r*','LineWidth',1)
ylabel('$y_{21}$ (volts)','FontSize',14,'Interpreter','Latex');
legend('Calculated h(t)','Measured h(t)')
grid on
xlabel('second','FontSize',14)
set(gca,'FontSize',14)
axis([0 2 -0.1 0.1])
hold off

subplot(414)
plot(t(42:N),h22d,'b*','LineWidth',1)
hold on
plot(t(42:N),y22(42:N),'r*','LineWidth',1)
ylabel('$y_{22}$ (volts)','FontSize',14,'Interpreter','Latex');
legend('Calculated h(t)','Measured h(t)')
grid on
xlabel('second','FontSize',14)
set(gca,'FontSize',14)
axis([0 2 -0.1 0.1])
hold off

%% Frequency Response / Ques #3

omega=linspace(0,20*2*pi,100);

for k=1:length(omega)
% For ns=6
temp1=C{1}*inv(exp(i*ts*omega(k))*eye(ns(1))-A{1})*B{1};    
freqrepa11(k)=temp1(1,1);
freqrepa12(k)=temp1(1,2);
freqrepa21(k)=temp1(2,1);
freqrepa22(k)=temp1(2,2);

% For ns=7
temp2=C{2}*inv(exp(i*ts*omega(k))*eye(ns(2))-A{2})*B{2};
freqrepb11(k)=temp2(1,1);
freqrepb12(k)=temp2(1,2);
freqrepb21(k)=temp2(2,1);
freqrepb22(k)=temp2(2,2);

% For ns=10
temp3=C{3}*inv(exp(i*ts*omega(k))*eye(ns(3))-A{3})*B{3};
freqrepc11(k)=temp3(1,1);
freqrepc12(k)=temp3(1,2);
freqrepc21(k)=temp3(2,1);
freqrepc22(k)=temp3(2,2);

% For ns=20
temp4=C{4}*inv(exp(i*ts*omega(k))*eye(ns(4))-A{4})*B{4};
freqrepd11(k)=temp4(1,1);
freqrepd12(k)=temp4(1,2);
freqrepd21(k)=temp4(2,1);
freqrepd22(k)=temp4(2,2);
end

figure(8)
loglog(omega/(2*pi),abs(freqrepa11),'LineWidth',1.3)
hold on
loglog(omega/(2*pi),abs(freqrepb11),'LineWidth',1.3)
loglog(omega/(2*pi),abs(freqrepc11),'LineWidth',1.3)
loglog(omega/(2*pi),abs(freqrepd11),'LineWidth',1.3)
ylabel('Magnitude [Volts]','FontSize',12.5,'Interpreter','Latex');
legend('$n_s=6$','$n_s=7$','$n_s=10$','$n_s=20$','FontSize',12.5,'Interpreter','Latex')
grid on
xlabel('Frequency [Hz]','FontSize',12.5,'Interpreter','Latex')
set(gca,'FontSize',12.5)
axis([0 20 1e-2 1])
hold off

figure(9)
semilogx(omega/(2*pi),rad2deg(angle(freqrepa11)),'LineWidth',1.3)
hold on
semilogx(omega/(2*pi),rad2deg(angle(freqrepb11)),'LineWidth',1.3)
semilogx(omega/(2*pi),rad2deg(angle(freqrepc11)),'LineWidth',1.3)
semilogx(omega/(2*pi),rad2deg(angle(freqrepd11)),'LineWidth',1.3)
ylabel('Phase [$^\circ$]','FontSize',12.5,'Interpreter','Latex');
legend('$n_s=6$','$n_s=7$','$n_s=10$','$n_s=20$','FontSize',12.5,'Interpreter','Latex')
grid on
xlabel('Frequency[Hz]','FontSize',12.5,'Interpreter','Latex')
set(gca,'FontSize',12.5)
axis([0 20 -200 200])
hold off

figure(10)
loglog(omega/(2*pi),abs(freqrepa12),'LineWidth',1.3)
hold on
loglog(omega/(2*pi),abs(freqrepb12),'LineWidth',1.3)
loglog(omega/(2*pi),abs(freqrepc12),'LineWidth',1.3)
loglog(omega/(2*pi),abs(freqrepd12),'LineWidth',1.3)
ylabel('Magnitude [Volts]','FontSize',12.5,'Interpreter','Latex');
legend('$n_s=6$','$n_s=7$','$n_s=10$','$n_s=20$','FontSize',12.5,'Interpreter','Latex')
grid on
xlabel('Frequency [Hz]','FontSize',12.5,'Interpreter','Latex')
set(gca,'FontSize',12.5)
axis([0 20 1e-2 1])
hold off

figure(11)
semilogx(omega/(2*pi),rad2deg(angle(freqrepa12)),'LineWidth',1.3)
hold on
semilogx(omega/(2*pi),rad2deg(angle(freqrepb12)),'LineWidth',1.3)
semilogx(omega/(2*pi),rad2deg(angle(freqrepc12)),'LineWidth',1.3)
semilogx(omega/(2*pi),rad2deg(angle(freqrepd12)),'LineWidth',1.3)
ylabel('Phase [$^\circ$]','FontSize',12.5,'Interpreter','Latex');
legend('$n_s=6$','$n_s=7$','$n_s=10$','$n_s=20$','FontSize',12.5,'Interpreter','Latex')
grid on
xlabel('Frequency[Hz]','FontSize',12.5,'Interpreter','Latex')
set(gca,'FontSize',12.5)
axis([0 20 -200 200])
hold off

figure(12)
loglog(omega/(2*pi),abs(freqrepa21),'LineWidth',1.3)
hold on
loglog(omega/(2*pi),abs(freqrepb21),'LineWidth',1.3)
loglog(omega/(2*pi),abs(freqrepc21),'LineWidth',1.3)
loglog(omega/(2*pi),abs(freqrepd21),'LineWidth',1.3)
ylabel('Magnitude [Volts]','FontSize',12.5,'Interpreter','Latex');
legend('$n_s=6$','$n_s=7$','$n_s=10$','$n_s=20$','FontSize',12.5,'Interpreter','Latex')
grid on
xlabel('Frequency [Hz]','FontSize',12.5,'Interpreter','Latex')
set(gca,'FontSize',12.5)
axis([0 20 1e-2 1])
hold off

figure(13)
semilogx(omega/(2*pi),rad2deg(angle(freqrepa21)),'LineWidth',1.3)
hold on
semilogx(omega/(2*pi),rad2deg(angle(freqrepb21)),'LineWidth',1.3)
semilogx(omega/(2*pi),rad2deg(angle(freqrepc21)),'LineWidth',1.3)
semilogx(omega/(2*pi),rad2deg(angle(freqrepd21)),'LineWidth',1.3)
ylabel('Phase [$^\circ$]','FontSize',12.5,'Interpreter','Latex');
legend('$n_s=6$','$n_s=7$','$n_s=10$','$n_s=20$','FontSize',12.5,'Interpreter','Latex')
grid on
xlabel('Frequency[Hz]','FontSize',12.5,'Interpreter','Latex')
set(gca,'FontSize',12.5)
axis([0 20 -200 200])
hold off

figure(14)
loglog(omega/(2*pi),abs(freqrepa22),'LineWidth',1.3)
hold on
loglog(omega/(2*pi),abs(freqrepb22),'LineWidth',1.3)
loglog(omega/(2*pi),abs(freqrepc22),'LineWidth',1.3)
loglog(omega/(2*pi),abs(freqrepd22),'LineWidth',1.3)
ylabel('Magnitude [Volts]','FontSize',12.5,'Interpreter','Latex');
legend('$n_s=6$','$n_s=7$','$n_s=10$','$n_s=20$','FontSize',12.5,'Interpreter','Latex')
grid on
xlabel('Frequency [Hz]','FontSize',12.5,'Interpreter','Latex')
set(gca,'FontSize',12.5)
axis([0 20 1e-2 1])
hold off

figure(15)
semilogx(omega/(2*pi),rad2deg(angle(freqrepa22)),'LineWidth',1.3)
hold on
semilogx(omega/(2*pi),rad2deg(angle(freqrepb22)),'LineWidth',1.3)
semilogx(omega/(2*pi),rad2deg(angle(freqrepc22)),'LineWidth',1.3)
semilogx(omega/(2*pi),rad2deg(angle(freqrepd22)),'LineWidth',1.3)
ylabel('Phase [$^\circ$]','FontSize',12.5,'Interpreter','Latex');
legend('$n_s=6$','$n_s=7$','$n_s=10$','$n_s=20$','FontSize',12.5,'Interpreter','Latex')
grid on
xlabel('Frequency[Hz]','FontSize',12.5,'Interpreter','Latex')
set(gca,'FontSize',12.5)
axis([0 20 -200 200])
hold off


%% Q4:Fast Fourier Transform
y11f = fft(y11(mi:N))./fft(u1(mi:N));
y12f = fft(y11(mi:N))./fft(u2(mi:N));
y21f = fft(y21(mi:N))./fft(u1(mi:N));
y22f = fft(y22(mi:N))./fft(u2(mi:N));
N = length(y11f);
om = [0:(361)-1]/(ts*(361));  %%%% frequency vector in hertz

figure(16)
loglog(omega/(2*pi),abs(freqrepa11),'LineWidth',1.3)
hold on
loglog(omega/(2*pi),abs(freqrepb11),'LineWidth',1.3)
loglog(omega/(2*pi),abs(freqrepc11),'LineWidth',1.3)
loglog(omega/(2*pi),abs(freqrepd11),'LineWidth',1.3)
loglog(om,abs(y11f),'LineWidth',1.3)
ylabel('Magnitude [Volts]','FontSize',12.5,'Interpreter','Latex');
legend('$n_s=6$','$n_s=7$','$n_s=10$','$n_s=20$','Emperical','FontSize',12.5,'Interpreter','Latex')
grid on
xlabel('Frequency [Hz]','FontSize',12.5,'Interpreter','Latex')
set(gca,'FontSize',12.5)
axis([0 20 1e-2 1])
hold off

figure(17)
semilogx(omega/(2*pi),rad2deg(angle(freqrepa11)),'LineWidth',1.3)
hold on
semilogx(omega/(2*pi),rad2deg(angle(freqrepb11)),'LineWidth',1.3)
semilogx(omega/(2*pi),rad2deg(angle(freqrepc11)),'LineWidth',1.3)
semilogx(omega/(2*pi),rad2deg(angle(freqrepd11)),'LineWidth',1.3)
semilogx(om,rad2deg(angle(y11f)),'LineWidth',1.3)
ylabel('Phase [$^\circ$]','FontSize',12.5,'Interpreter','Latex');
legend('$n_s=6$','$n_s=7$','$n_s=10$','$n_s=20$','Emperical','FontSize',12.5,'Interpreter','Latex')
grid on
xlabel('Frequency[Hz]','FontSize',12.5,'Interpreter','Latex')
set(gca,'FontSize',12.5)
axis([0 20 -200 200])
hold off

figure(18)
loglog(omega/(2*pi),abs(freqrepa12),'LineWidth',1.3)
hold on
loglog(omega/(2*pi),abs(freqrepb12),'LineWidth',1.3)
loglog(omega/(2*pi),abs(freqrepc12),'LineWidth',1.3)
loglog(omega/(2*pi),abs(freqrepd12),'LineWidth',1.3)
loglog(om,abs(y12f),'LineWidth',1.3)
ylabel('Magnitude [Volts]','FontSize',12.5,'Interpreter','Latex');
legend('$n_s=6$','$n_s=7$','$n_s=10$','$n_s=20$','Emperical','FontSize',12.5,'Interpreter','Latex')
grid on
xlabel('Frequency [Hz]','FontSize',12.5,'Interpreter','Latex')
set(gca,'FontSize',12.5)
axis([0 20 1e-2 1])
hold off

figure(19)
semilogx(omega/(2*pi),rad2deg(angle(freqrepa12)),'LineWidth',1.3)
hold on
semilogx(omega/(2*pi),rad2deg(angle(freqrepb12)),'LineWidth',1.3)
semilogx(omega/(2*pi),rad2deg(angle(freqrepc12)),'LineWidth',1.3)
semilogx(omega/(2*pi),rad2deg(angle(freqrepd12)),'LineWidth',1.3)
semilogx(om,rad2deg(angle(y12f)),'LineWidth',1.3)
ylabel('Phase [$^\circ$]','FontSize',12.5,'Interpreter','Latex');
legend('$n_s=6$','$n_s=7$','$n_s=10$','$n_s=20$','Emperical','FontSize',12.5,'Interpreter','Latex')
grid on
xlabel('Frequency[Hz]','FontSize',12.5,'Interpreter','Latex')
set(gca,'FontSize',12.5)
axis([0 20 -200 200])
hold off

figure(20)
loglog(omega/(2*pi),abs(freqrepa21),'LineWidth',1.3)
hold on
loglog(omega/(2*pi),abs(freqrepb21),'LineWidth',1.3)
loglog(omega/(2*pi),abs(freqrepc21),'LineWidth',1.3)
loglog(omega/(2*pi),abs(freqrepd21),'LineWidth',1.3)
loglog(om,abs(y21f),'LineWidth',1.3)
ylabel('Magnitude [Volts]','FontSize',12.5,'Interpreter','Latex');
legend('$n_s=6$','$n_s=7$','$n_s=10$','$n_s=20$','Emperical','FontSize',12.5,'Interpreter','Latex')
grid on
xlabel('Frequency [Hz]','FontSize',12.5,'Interpreter','Latex')
set(gca,'FontSize',12.5)
axis([0 20 1e-2 1])
hold off

figure(21)
semilogx(omega/(2*pi),rad2deg(angle(freqrepa21)),'LineWidth',1.3)
hold on
semilogx(omega/(2*pi),rad2deg(angle(freqrepb21)),'LineWidth',1.3)
semilogx(omega/(2*pi),rad2deg(angle(freqrepc21)),'LineWidth',1.3)
semilogx(omega/(2*pi),rad2deg(angle(freqrepd21)),'LineWidth',1.3)
semilogx(om,rad2deg(angle(y21f)),'LineWidth',1.3)
ylabel('Phase [$^\circ$]','FontSize',12.5,'Interpreter','Latex');
legend('$n_s=6$','$n_s=7$','$n_s=10$','$n_s=20$','Emperical','FontSize',12.5,'Interpreter','Latex')
grid on
xlabel('Frequency[Hz]','FontSize',12.5,'Interpreter','Latex')
set(gca,'FontSize',12.5)
axis([0 20 -200 200])
hold off

figure(22)
loglog(omega/(2*pi),abs(freqrepa22),'LineWidth',1.3)
hold on
loglog(omega/(2*pi),abs(freqrepb22),'LineWidth',1.3)
loglog(omega/(2*pi),abs(freqrepc22),'LineWidth',1.3)
loglog(omega/(2*pi),abs(freqrepd22),'LineWidth',1.3)
loglog(om,abs(y22f),'LineWidth',1.3)
ylabel('Magnitude [Volts]','FontSize',12.5,'Interpreter','Latex');
legend('$n_s=6$','$n_s=7$','$n_s=10$','$n_s=20$','Emperical','FontSize',12.5,'Interpreter','Latex')
grid on
xlabel('Frequency [Hz]','FontSize',12.5,'Interpreter','Latex')
set(gca,'FontSize',12.5)
axis([0 20 1e-2 1])
hold off

figure(23)
semilogx(omega/(2*pi),rad2deg(angle(freqrepa22)),'LineWidth',1.3)
hold on
semilogx(omega/(2*pi),rad2deg(angle(freqrepb22)),'LineWidth',1.3)
semilogx(omega/(2*pi),rad2deg(angle(freqrepc22)),'LineWidth',1.3)
semilogx(omega/(2*pi),rad2deg(angle(freqrepd22)),'LineWidth',1.3)
semilogx(om,rad2deg(angle(y22f)),'LineWidth',1.3)
ylabel('Phase [$^\circ$]','FontSize',12.5,'Interpreter','Latex');
legend('$n_s=6$','$n_s=7$','$n_s=10$','$n_s=20$','Emperical','FontSize',12.5,'Interpreter','Latex')
grid on
xlabel('Frequency[Hz]','FontSize',12.5,'Interpreter','Latex')
set(gca,'FontSize',12.5)
axis([0 20 -200 200])
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Task 2 / Ques #1

sys = ss(A{2},B{2},C{2},0,ts);  % Defining discrete State Space system
% Calculating the trans zero using build in function
transzero=tzero(sys) 

% Solving Transmission zero using generalized eigenvalue method
TZa=[A{2} B{2};
    -C{2} zeros(2)];
TZb=[eye(7) zeros(7,2);
    zeros(2,9)];

zero=eig(TZa,TZb)  % Transmission Zero
eigA=eig(A{2});     % Eigenvalues of A with ns=7

%% Ques #2

circ = exp(1i*[0:360]*pi/180);  % Crete circle for compairison in plot

figure(24)
plot(real(eigA),imag(eigA),'x','LineWidth',1.3, 'MarkerSize',8)
hold on
grid on
plot(real(zero),imag(zero),'o','LineWidth',1.3, 'MarkerSize',8)
plot(real(circ),imag(circ),'r--','LineWidth',1.1)
axis square
legend('Eigenvalue','Zeros','$R_{circle}=1$','FontSize',12.5,'Interpreter','Latex')
xlabel('Real Part','FontSize',12.5,'Interpreter','Latex')
ylabel('Imaginary Part','FontSize',12.5,'Interpreter','Latex')
axis([-3 3 -3 3])
hold off

%% Ques #3

lamd_c=log(eigA)/ts
damped=lamd_c(1:4)
nat_freq=abs(damped)

%% Ques #4

B1=B{2}(:,1);   % Closing the 2nd input
B2=B{2}(:,2);   % Closing the 1st input
C1=C{2}(1,:);   % Closing the 2nd output
C2=C{2}(2,:);   % Closing the 1st output

% Case 1 (Input 1 and Output 1 are on)
TZac1=[A{2} B1; -C1 0];
TZbc1=[eye(7) zeros(7,1); zeros(1,8)];
zeroc1=eig(TZac1,TZbc1); % Generalized Eigenvalue Problem

% Case 2 (Input 1 and Output 2 are on)
TZac2=[A{2} B1; -C2 0];   
TZbc2=[eye(7) zeros(7,1); zeros(1,8)];
zeroc2=eig(TZac2,TZbc2); % Generalized Eigenvalue Problem

% Case 3 (Input 2 and Output 1 are on)
TZac3=[A{2} B2; -C1 0];
TZbc3=[eye(7) zeros(7,1); zeros(1,8)];
zeroc3=eig(TZac3,TZbc3); % Generalized Eigenvalue Problem

% Case 4 (Input 2 and Output 2 are on)
TZac4=[A{2} B2; -C2 0];
TZbc4=[eye(7) zeros(7,1); zeros(1,8)];
zeroc4=eig(TZac4,TZbc4); % Generalized Eigenvalue Problem

figure(25)
plot(real(eigA),imag(eigA),'x','LineWidth',1.3, 'MarkerSize',8)
hold on
grid on
plot(real(zeroc1),imag(zeroc1),'o','LineWidth',1.3, 'MarkerSize',8)
plot(real(circ),imag(circ),'r--','LineWidth',1.1)
axis square
legend('Eigenvalue','Zeros','$R_{circle}=1$','FontSize',12.5,'Interpreter','Latex')
xlabel('Real Part','FontSize',12.5,'Interpreter','Latex')
ylabel('Imaginary Part','FontSize',12.5,'Interpreter','Latex')
axis([-3 3 -3 3])
hold off

figure(26)
plot(real(eigA),imag(eigA),'x','LineWidth',1.3, 'MarkerSize',8)
hold on
grid on
plot(real(zeroc2),imag(zeroc2),'o','LineWidth',1.3, 'MarkerSize',8)
plot(real(circ),imag(circ),'r--','LineWidth',1.1)
axis square
legend('Eigenvalue','Zeros','$R_{circle}=1$','FontSize',12.5,'Interpreter','Latex')
xlabel('Real Part','FontSize',12.5,'Interpreter','Latex')
ylabel('Imaginary Part','FontSize',12.5,'Interpreter','Latex')
axis([-3 3 -3 3])
hold off

figure(27)
plot(real(eigA),imag(eigA),'x','LineWidth',1.3, 'MarkerSize',8)
hold on
grid on
plot(real(zeroc3),imag(zeroc3),'o','LineWidth',1.3, 'MarkerSize',8)
plot(real(circ),imag(circ),'r--','LineWidth',1.1)
axis square
legend('Eigenvalue','Zeros','$R_{circle}=1$','FontSize',12.5,'Interpreter','Latex')
xlabel('Real Part','FontSize',12.5,'Interpreter','Latex')
ylabel('Imaginary Part','FontSize',12.5,'Interpreter','Latex')
axis([-3 3 -3 3])
hold off

figure(28)
plot(real(eigA),imag(eigA),'x','LineWidth',1.3, 'MarkerSize',8)
hold on
grid on
plot(real(zeroc4),imag(zeroc4),'o','LineWidth',1.3, 'MarkerSize',8)
plot(real(circ),imag(circ),'r--','LineWidth',1.1)
axis square
legend('Eigenvalue','Zeros','$R_{circle}=1$','FontSize',12.5,'Interpreter','Latex')
xlabel('Real Part','FontSize',12.5,'Interpreter','Latex')
ylabel('Imaginary Part','FontSize',12.5,'Interpreter','Latex')
axis([-3 3 -3 3])
hold off

% Hankel Matrices for each channel
for q=1:100
for r=1:100 
    H11(q,r)=y11(r+q+mi-1);
    H12(q,r)=y12(r+q+mi-1);
    H21(q,r)=y21(r+q+mi-1);
    H22(q,r)=y22(r+q+mi-1);
end
end

svdH11=svd(H11);
svdH12=svd(H12);
svdH21=svd(H21);
svdH22=svd(H22);

figure(29)
subplot(221)
semilogy(svdx(1:100),svdH11,'o', 'MarkerFaceColor', 'b')
grid on
axis([0 40 1e-3 1])
xlabel('Singular Value Index')
ylabel('Hankel Singular Value')
title('H_{y_{11}} Singular Values for the First 40')

subplot(222)
semilogy(svdx(1:100),svdH12,'o', 'MarkerFaceColor', 'b')
grid on
axis([0 40 1e-3 1])
xlabel('Singular Value Index')
ylabel('Hankel Singular Value')
title('H_{y_{12}} Singular Values for the First 40')

subplot(223)
semilogy(svdx(1:100),svdH21,'o', 'MarkerFaceColor', 'b')
grid on
axis([0 40 1e-3 1])
xlabel('Singular Value Index')
ylabel('Hankel Singular Value')
title('H_{y_{21}} Singular Values for the First 40')

subplot(224)
semilogy(svdx(1:100),svdH22,'o', 'MarkerFaceColor', 'b')
grid on
axis([0 40 1e-3 1])
xlabel('Singular Value Index')
ylabel('Hankel Singular Value')
title('H_{y_{22}} Singular Values for the First 40')

%% Ques #5 for ns=8

eigA=eig(A{5});  % Eigval of A with ns=8
B1=B{5}(:,1);    % Closing the 2nd input
B2=B{5}(:,2);    % Closing the 1st input
C1=C{5}(1,:);    % Closing the 2nd output
C2=C{5}(2,:);    % Closing the 1st output

% Case 1 (Input 1 and Output 1 are on)
TZac1=[A{5} B1; -C1 0];
TZbc1=[eye(8) zeros(8,1); zeros(1,9)];
zeroc1=eig(TZac1,TZbc1); % Solve for zeros using gen eigval prob

% Case 2 (Input 1 and Output 2 are on)
TZac2=[A{5} B1;  -C2 0];
TZbc2=[eye(8) zeros(8,1); zeros(1,9)];
zeroc2=eig(TZac2,TZbc2); % Solve for zeros using gen eigval prob

% Case 3 (Input 2 and Output 1 are on)
TZac3=[A{5} B2; -C1 0];
TZbc3=[eye(8) zeros(8,1); zeros(1,9)];
zeroc3=eig(TZac3,TZbc3); % Solve for zeros using gen eigval prob

% Case 4 (Input 2 and Output 2 are on)
TZac4=[A{5} B2; -C2 0];    
TZbc4=[eye(8) zeros(8,1); zeros(1,9)];
zeroc4=eig(TZac4,TZbc4); % Solve for zeros using gen eigval prob

figure(30)
plot(real(eigA),imag(eigA),'x','LineWidth',1.3, 'MarkerSize',8)
hold on
grid on
plot(real(zeroc1),imag(zeroc1),'o','LineWidth',1.3, 'MarkerSize',8)
plot(real(circ),imag(circ),'r--','LineWidth',1.1)
axis square
legend('Eigenvalue','Zeros','$R_{circle}=1$','FontSize',12.5,'Interpreter','Latex')
xlabel('Real Part','FontSize',12.5,'Interpreter','Latex')
ylabel('Imaginary Part','FontSize',12.5,'Interpreter','Latex')
axis([-3 3 -3 3])
hold off

figure(31)
plot(real(eigA),imag(eigA),'x','LineWidth',1.3, 'MarkerSize',8)
hold on
grid on
plot(real(zeroc2),imag(zeroc2),'o','LineWidth',1.3, 'MarkerSize',8)
plot(real(circ),imag(circ),'r--','LineWidth',1.1)
axis square
legend('Eigenvalue','Zeros','$R_{circle}=1$','FontSize',12.5,'Interpreter','Latex')
xlabel('Real Part','FontSize',12.5,'Interpreter','Latex')
ylabel('Imaginary Part','FontSize',12.5,'Interpreter','Latex')
axis([-3 3 -3 3])
hold off

figure(32)
plot(real(eigA),imag(eigA),'x','LineWidth',1.3, 'MarkerSize',8)
hold on
grid on
plot(real(zeroc3),imag(zeroc3),'o','LineWidth',1.3, 'MarkerSize',8)
plot(real(circ),imag(circ),'r--','LineWidth',1.1)
axis square
legend('Eigenvalue','Zeros','$R_{circle}=1$','FontSize',12.5,'Interpreter','Latex')
xlabel('Real Part','FontSize',12.5,'Interpreter','Latex')
ylabel('Imaginary Part','FontSize',12.5,'Interpreter','Latex')
axis([-3 3 -3 3])
hold off

figure(33)
plot(real(eigA),imag(eigA),'x','LineWidth',1.3, 'MarkerSize',8)
hold on
grid on
plot(real(zeroc4),imag(zeroc4),'o','LineWidth',1.3, 'MarkerSize',8)
plot(real(circ),imag(circ),'r--','LineWidth',1.1)
axis square
legend('Eigenvalue','Zeros','$R_{circle}=1$','FontSize',12.5,'Interpreter','Latex')
xlabel('Real Part','FontSize',12.5,'Interpreter','Latex')
ylabel('Imaginary Part','FontSize',12.5,'Interpreter','Latex')
axis([-3 3 -3 3])
hold off