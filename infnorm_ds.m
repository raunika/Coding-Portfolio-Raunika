function [gam omeg]=infnorm_ds(A,B,C,D,up_lim,lo_lim,tol)
% 
% 
% The function output is gamma/ infinity norm with the frequency occured
% up_lim is the upper bound of gamma
% lo_lim is the lower bound of gamma
% tol is the tolerance or accuracy of gamma required

% Convert Discrete time SS to Continous time SS
Ac=-inv(eye(length(A))+A)*(eye(length(A))-A);
Bc=sqrt(2)*inv(eye(length(A))+A)*B;
Cc=sqrt(2)*C*inv(eye(length(A))+A);
Dc=D-C*inv(eye(length(A))+A)*B;

% Using the Function developed for continous ss to calculate infinity norm
[gam omeg]=infnorm_cs(Ac,Bc,Cc,Dc,up_lim,lo_lim,tol);




end