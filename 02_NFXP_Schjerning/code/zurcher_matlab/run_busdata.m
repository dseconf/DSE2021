% run_busdata: Estimates Rust's engine repplacement model using the Busdata from Rust(Ecta, 1987) 
% clear

global ev0 V0;  
ev0=0; V0=0;

%% Parameters for solution algorithm (used insolve.m)
% ap=solve.setup;

% Swiches:
twostep=0;      	% 1: Estimate using two-step pmle, 0: Use full mle
bustypes=4;		  	% Select bus types smaller than this number, i.e. 1,2,..bustypes

% Read default parameters in to struct mp
mpopt.n=175;
mpopt.integrated=0;
mp=zurcher.setup(mpopt);

% Starting values for parameters to be estimated
mp.RC=0;     		% Replacement cost
mp.c=0;				% Cost parameter

data = zurcher.readbusdata(mp, bustypes);

[nfxp_results, pnames, theta_hat, Avar]=nfxp.estim(data, mp, twostep);

% ************************************
% Print output
% ************************************
fprintf('Structural Estimation using busdata from Rust(1987)\n');
fprintf('Beta           = %10.5f \n',mp.beta);
fprintf('n              = %10.5f \n',mp.n);
fprintf('Sample size    = %10.5f \n',numel(data.d));
fprintf('\n'); 
fprintf('\n'); 

mphat = output.estimates(nfxp_results, pnames, theta_hat, Avar);
fprintf('log-likelihood    = %10.5f \n',nfxp_results.llval);
fprintf('runtime (seconds) = %10.5f \n',nfxp_results.cputime);
