global ev0;  ev0 = 0;

%% Parameters for solution algorithm (used in solve.m)
ap=solve.setup;

% ************************************
% Model solution used for DGP
% ************************************
% Swiches:
twostep=1;      % 1: Estimate using two-step pmle (FULL MLE is not implemented for MPEC)

nMC=10; 	  	% number of MC samples
N=50;			  	% Number of busses to simulate 
T=119;			% Number of time periods to simulate 

% Spaces
mp.n=175;					% Number of gridpoints		
mp.max=450;				% Max of mileage
mp.grid= (0:mp.n-1)';  	% Grid over mileage

% Structural parameters
mp.p=[0.0937 0.4475 0.4459 0.0127]';   	% Transition probabiliuties
mp.RC=11.7257;     											% Replacement cost
mp.c=2.45569;														% Cost parameter
mp.beta=0.9999;													% Discount factor

% Read default parameters in to struct mp 
mpopt.integrated=0;
mp=zurcher.setup(mpopt);	


% Cost function
cost0=0.001*mp.c*mp.grid;				

% Transition matrix for mileage
P0 = zurcher.statetransition(mp);	

bellman= @(ev) zurcher.bellman_ev(ev, mp, P0);
[~, pk0]=solve.poly(bellman, ev0, ap, mp.beta);	

% ************************************
% START MONTE CARLO HERE
% ************************************

fprintf('Begin Monte Carlo, with n=%d replications\n', nMC);
rand('seed',300);

mpec_results=struct;
nfxp_results=struct;

%main loop
for i_mc = 1:nMC;
	% ************************************
	% STEP 1: SIMULATE DATA 
	% ************************************

	timetosimulate=tic;
	data = zurcher.simdata(N, T, mp, P0, pk0);
	timetosimulate=toc(timetosimulate);
	fprintf('i_mc=%d, Time to simulate data : %1.5g seconds\n', i_mc, timetosimulate);

  	% ************************************
	% STEP 2a: ESTIMATE parameters using NFXP
	% ************************************

	[nfxp_result_i, pnames, theta_hat, Avar]=nfxp.estim(data, mp, twostep);
	nfxp_results=output.savemc(nfxp_result_i, nfxp_results, i_mc);

	% ************************************
 	% STEP 2b: ESTIMATE parameters using MPEC
 	% ************************************

	[mpec_result_i, pnames, theta_hat]=mpec.estim(data, mp, twostep);
	mpec_results=output.savemc(mpec_result_i, mpec_results, i_mc);

end  % End Monte Carlo

nfxp_results.title = 'Nestec Fixed Point Algorithm, NFXP';
mpec_results.title = 'Mathematical Programming with Equilibrium Constraints (MPEC)';

results = {nfxp_results,mpec_results};

output.table_mc({'RC', 'c'}, mp, results); 
output.table_np(results); 
