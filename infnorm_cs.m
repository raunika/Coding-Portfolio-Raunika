function [gam omeg]=infnorm_cs(A,B,C,D,up_lim,lo_lim,tol)

% The function output is gamma/ infinity norm with the frequency occured
% up_lim is the upper bound of gamma
% lo_lim is the lower bound of gamma
% tol is the tolerance or accuracy of gamma required

%% Parameters
n=length(D);        % number of input/output
ts=1/40;            % Only applicable to this experiment

%% Check and update upper and lower bound

while (up_lim-lo_lim)/2 > tol % Running calculation if within tolerence

    
gam=(lo_lim+up_lim)/2;      % Select gam as middle value and start checking
D_gam=gam^2*eye(n)-D'*D;    % D_gamma

% Construct A_clp matrix
A_clp=[A+B*inv(D_gam)*D'*C  -B*inv(D_gam)*B';
       C'*C+C'*D*inv(D_gam)*D'*C  -A'-C'*D*inv(D_gam)*B'];
   
eval=eig(A_clp);        % Eigenvalues of A_clp
eval_real=real(eval);   % Real Part of Eigenvalues A_clp
eval_imag=imag(eval);   % Imag Part of Eigenvalues A_clp
nu=max(eval_imag);      % For the frequency calculation

for k=1:length(eval) % checking if purely imag eval exists
    if (abs(eval_real(k))<1e-8 && abs(eval_imag(k))>0) 
        lo_lim=gam; % gam is the lower bound
    end
end 

% gam is the upper bound if purely imag eval doesnt exists
if lo_lim~=gam
    up_lim=gam;
end

end
  
%% Calculate Omega and finalize gamma

exp_jomegts=(1+j*nu)/(1-j*nu);
omeg=real(log(exp_jomegts)/(j*ts))/(2*pi); % unit Hz
gam=(lo_lim+up_lim)/2; 

end