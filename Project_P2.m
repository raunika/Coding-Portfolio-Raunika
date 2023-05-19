clear all;
close all;
clc

%% Loading Data
load u_rand.mat
y1 = u_rand.Y(3).Data;
y2 = u_rand.Y(4).Data;
u1 = u_rand.Y(1).Data;
u2 = u_rand.Y(2).Data;
ts = 1/40;
N = length(y1);
t = [0:N-1]*ts - 1;
u=[u1; u2];

load u1_impulse.mat

y11 = u1_impulse.Y(3).Data;
y21 = u1_impulse.Y(4).Data;
u11 = u1_impulse.Y(1).Data;  %%% note that the pulse magnitude is 5

[m,mi] = max(u1>0);  %%% find index where pulse occurs

load u2_impulse.mat

y12 = u2_impulse.Y(3).Data;
y22 = u2_impulse.Y(4).Data;
u22 = u2_impulse.Y(2).Data;

%%% remove any offsets in output data using data prior to pulse application
y11 = y11 - mean(y11([1:mi-1]));
y12 = y12 - mean(y12([1:mi-1]));
y21 = y21 - mean(y21([1:mi-1]));
y22 = y22 - mean(y22([1:mi-1]));

%%% rescale IO data so that impulse input has magnitude 1
y11 = y11/max(u11);
y12 = y12/max(u22);
y21 = y21/max(u11);
y22 = y22/max(u22);
u11 = u11/max(u11);
u22 = u22/max(u22);

NN = length(y11);  %%%% length of data sets
tt = [0:NN-1]*ts - 1;

%% Ques #1

mean(u1)    % Verify mean of input approx zero
mean(u2)    % Verify mean of input approx zero

%% Ques #2

% Setting up variable for auto-correlation
Ruu11=zeros(1,401);
Ruu21=zeros(1,401);
Ruu12=zeros(1,401);
Ruu22=zeros(1,401);

% Calculate for Ruu
for k=-200:200 % from t=-5 to t=5
    
    for q=-11800:11800 
        
        temp=[u1(k+q+12001) ; u2(k+q+12001)]*[u1(q+12001) u2(q+12001)];
        Ruu11(k+201)=Ruu11(k+201)+temp(1,1); 
        Ruu21(k+201)=Ruu21(k+201)+temp(2,1);
        Ruu12(k+201)=Ruu12(k+201)+temp(1,2);
        Ruu22(k+201)=Ruu22(k+201)+temp(2,2);

    end
    Ruu11(k+201)=(1/(2*11800+1))*Ruu11(k+201); 
    Ruu21(k+201)=(1/(2*11800+1))*Ruu21(k+201);
    Ruu12(k+201)=(1/(2*11800+1))*Ruu12(k+201);
    Ruu22(k+201)=(1/(2*11800+1))*Ruu22(k+201);
end

% Create array for time
tau_uu=linspace(-5,5,401);

figure(1)
plot(tau_uu,Ruu11,'o','MarkerFaceColor', 'r')
grid on
axis([-5 5 -0.5 4.5])
xlabel('Time [s]','FontSize',12.5,'Interpreter','Latex')
% ylabel('Magnitude','FontSize',12.5,'Interpreter','Latex')
title('$R_{uu_{11}}$ Plot for [-5 5]','FontSize',12.5,'Interpreter','Latex')

figure(2)
plot(tau_uu,Ruu12,'o','MarkerFaceColor', 'r')
grid on
axis([-5 5 -0.5 4.5])
xlabel('Time [s]','FontSize',12.5,'Interpreter','Latex')
% ylabel('Magnitude','FontSize',12.5,'Interpreter','Latex')
title('$R_{uu_{12}}$ Plot for [-5 5]','FontSize',12.5,'Interpreter','Latex')

figure(3)
plot(tau_uu,Ruu21,'o','MarkerFaceColor', 'r')
grid on
axis([-5 5 -0.5 4.5])
xlabel('Time [s]','FontSize',12.5,'Interpreter','Latex')
% ylabel('Magnitude','FontSize',12.5,'Interpreter','Latex')
title('$R_{uu_{21}}$ Plot for [-5 5]','FontSize',12.5,'Interpreter','Latex')

figure(4)
plot(tau_uu,Ruu22,'o','MarkerFaceColor', 'r')
grid on
axis([-5 5 -0.5 4.5])
xlabel('Time [s]','FontSize',12.5,'Interpreter','Latex')
% ylabel('Magnitude','FontSize',12.5,'Interpreter','Latex')
title('$R_{uu_{22}}$ Plot for [-5 5]','FontSize',12.5,'Interpreter','Latex')

%% Ques #3
varu1=var(u1);
varu2=var(u2);

%% Ques #4

% Defining var for cross corelation
Ryu11=zeros(1,89);
Ryu21=zeros(1,89);
Ryu12=zeros(1,89);
Ryu22=zeros(1,89);

for k=-8:80 % from t=-0.2 to t=2
    
    for q=-11920:11920
        
        temp=[y1(k+q+12001) ; y2(k+q+12001)]*[u1(q+12001) u2(q+12001)];
        Ryu11(k+9)=Ryu11(k+9)+temp(1,1); 
        Ryu21(k+9)=Ryu21(k+9)+temp(2,1);
        Ryu12(k+9)=Ryu12(k+9)+temp(1,2);
        Ryu22(k+9)=Ryu22(k+9)+temp(2,2);

    end
    % Normalized by number of data and variance of input
    Ryu11(k+9)=(1/(2*11920+1))*Ryu11(k+9)*(1/varu1);  
    Ryu21(k+9)=(1/(2*11920+1))*Ryu21(k+9)*(1/varu1);
    Ryu12(k+9)=(1/(2*11920+1))*Ryu12(k+9)*(1/varu2);
    Ryu22(k+9)=(1/(2*11920+1))*Ryu22(k+9)*(1/varu2);
end

% Define time vector
tau_yu=linspace(-0.2,2,89);

figure(5)
plot(tau_yu,Ryu11,'o','MarkerFaceColor', 'b')
hold on
plot (tt,y11,'o','MarkerFaceColor', 'g')
grid on
axis([-0.2 2 -0.1 0.1])
xlabel('Time [s]','FontSize',12.5,'Interpreter','Latex')
% ylabel('Magnitude','FontSize',12.5,'Interpreter','Latex')
title('$R_{yu_{11}}$ Plot for [-0.2 2]','FontSize',12.5,'Interpreter','Latex')
legend('$R_{yu_{11}}$','$h_{11}$','FontSize',12.5,'Interpreter','Latex')
hold off

figure(6)
plot(tau_yu,Ryu12,'o','MarkerFaceColor', 'b')
hold on
plot (tt,y12,'o','MarkerFaceColor', 'g')
grid on
axis([-0.2 2 -0.1 0.1])
xlabel('Time [s]','FontSize',12.5,'Interpreter','Latex')
% ylabel('Magnitude','FontSize',12.5,'Interpreter','Latex')
title('$R_{yu_{12}}$ Plot for [-0.2 2]','FontSize',12.5,'Interpreter','Latex')
legend('$R_{yu_{12}}$','$h_{12}$','FontSize',12.5,'Interpreter','Latex')
hold off

figure(7)
plot(tau_yu,Ryu21,'o','MarkerFaceColor', 'b')
hold on
plot (tt,y21,'o','MarkerFaceColor', 'g')
grid on
axis([-0.2 2 -0.1 0.1])
xlabel('Time [s]','FontSize',12.5,'Interpreter','Latex')
% ylabel('Magnitude','FontSize',12.5,'Interpreter','Latex')
title('$R_{yu_{21}}$ Plot for [-0.2 2]','FontSize',12.5,'Interpreter','Latex')
legend('$R_{yu_{21}}$','$h_{21}$','FontSize',12.5,'Interpreter','Latex')
hold off

figure(8)
plot(tau_yu,Ryu22,'o','MarkerFaceColor', 'b')
hold on
plot (tt,y22,'o','MarkerFaceColor', 'g')
grid on
axis([-0.2 2 -0.1 0.1])
xlabel('Time [s]','FontSize',12.5,'Interpreter','Latex')
% ylabel('Magnitude','FontSize',12.5,'Interpreter','Latex')
title('$R_{yu_{22}}$ Plot for [-0.2 2]','FontSize',12.5,'Interpreter','Latex')
legend('$R_{yu_{22}}$','$h_{22}$','FontSize',12.5,'Interpreter','Latex')
hold off

