% run_busdata: Estimates Rust's engine repplacement model using the Busdata from Rust(Ecta, 1987) 
% Metods used in this script: NFXP full me
% clear

global ev0 V0;  
ev0=0; V0=0;

bustypes=4;		  	% Select bus types smaller than this number, i.e. 1,2,..bustypes

% methods to cosinder
method={'nfxp (mle)'}; % only Full MLE
% method={'nfxp (mle)', 'mpec (pmle)', 'nfxp (pmle)'}; % compare NFXP and MPEC

% Read default parameters in to struct mp
mpopt.n=175;
mpopt.integrated=0; % 1: use integrated value function, 0: use expected value function
mp=zurcher.setup(mpopt);

% Starting values for parameters to be estimated
mp.RC=0;     		% Replacement cost
mp.c=0;				% Cost parameter

data = zurcher.readbusdata(mp, bustypes);

for i=1:numel(method)
	switch method{i}
		case 'nfxp (mle)'
			% Full MLE using NFXP implementation 
			[results, pnames, theta_hat, Avar]=nfxp.estim(data, mp, 0);
		case 'nfxp (pmle)'
			% Two step partial MLE using NFXP
			[results, pnames, theta_hat, Avar]=nfxp.estim(data, mp, 1);
		case 'mpec (pmle)'
			% Two step partial MLE using MPEC
			[results, pnames, theta_hat, Avar]=mpec.estim(data, mp, 1);
		otherwise
			error('Method does not exist');
	end

	% ************************************
	% Print output
	% ************************************

	if i==1	
		fprintf('Structural Estimation using busdata from Rust(1987)\n');
		fprintf('Beta           = %10.5f \n',mp.beta);
		fprintf('n              = %10.5f \n',mp.n);
		fprintf('Sample size    = %10.5f \n',numel(data.d));
		fprintf('\n'); 
	end
	fprintf('\n'); 

	fprintf('\nMethod %s\n', method{i});


	mphat = output.estimates(results, pnames, theta_hat, Avar);
	fprintf('log-likelihood    = %10.5f \n',results.llval);
	fprintf('runtime (seconds) = %10.5f \n',results.cputime);
end







