% run_mc_nfxp: Monte Carlo experiement to asses the performance of NFXP on Rust's engine repplacement model

global ev0;  ev0=0;

%% Parameters for solution algorithm (used insolve.m)
ap=solve.setup;
ap.printfxp=0;

% Swiches:
twostep=0;      % 1: Estimate using two-step pmle, 0: Use full mle
nMC=10; 		% number of MC samples
N=50;			% Number of busses to simulate 
T=119;			% Number of time periods to simulate 

% Read default parameters in to struct mp 
mp=zurcher.setup;		

% Transition matrix for mileage
P0 = zurcher.statetransition(mp);	

% Model solution used for DGP
bellman= @(ev) zurcher.bellman_ev(ev, mp,  P0);
[~, pk0]=solve.poly(bellman, ev0, ap, mp.beta);	

fprintf('Begin Monte Carlo, with n=%d replications\n', nMC);

rand('seed',300);

nfxp_results=struct; 
for i_mc=1:nMC;
	% ************************************
	% SIMULATE DATA 
	% ************************************

	timetosimulate=tic;
	data = zurcher.simdata(N, T, mp, P0, pk0);
	timetosimulate=toc(timetosimulate);
	fprintf('i_mc=%d, Time to simulate data : %1.5g seconds\n', i_mc, timetosimulate);

	% ************************************
	% Estimate parameters and collect results
	% ************************************
	result_i=nfxp.estim(data, mp, twostep);
	nfxp_results=output.savemc(result_i, nfxp_results, i_mc);
end  % End Monte Carlo

nfxp_results.title = 'Monte Carlo results, NFXP';
output.table_mc({'RC', 'c'}, mp, {nfxp_results}); 
output.table_np({nfxp_results}); 

