clear;
close all
global ev0;  ev0=0;

do_nk=0; % 1 for Newton iterations, 0 for succesive approximations

%% Parameters for solution algorithm (used insolve.m)
ap=solve.setup;

% Read default parameters in to struct mp 
mp=zurcher.setup;
mp.integrated=0;

% Spaces
mp.n=90;				% Number of gridpoints		
mp.grid= (0:mp.n-1)';  	% Grid over mileage

% Structural parameters
mp.p=[0.0937 0.4475 0.4459 0.0127]';   	% Transition probabiliuties
mp.RC=11.7257;     						% Replacement cost
mp.c=2.45569;							% Cost parameter
mp.beta=0.9999;							% Discount factor

P0 = zurcher.statetransition(mp);	% Transition matrix for mileage

legends={};

if do_nk
	nsaiter=5;
else
	nsaiter=10000;
end

colorOrder = get(gca, 'ColorOrder');
ap.sa_min=nsaiter;
ap.sa_max=ap.sa_min;
ap.printfxp   = 0;            	% Print iteration info for fixed point algorithm if > 0
betavec=[0.95 0.99 0.999 0.9999];
for ibeta=1:numel(betavec)
	mp.beta=betavec(ibeta);
  legends{ibeta}=sprintf('beta=%g', betavec(ibeta)) ;

  % belman equation
	bellman= @(ev) zurcher.bellman(ev, mp, P0);
	[ev, pk0, dev, iterinfo(ibeta)]=solve.poly(bellman, ev0, ap, mp.beta);	


	% [mp.ev, pk0, F, iterinfo(ibeta)]=nfxp.solve(0, P0, cost0, mp, nfxp_options); 	% Solve the model
	tol=iterinfo(ibeta).sa.tol;

	if do_nk
		tol=[tol; iterinfo(ibeta).nk.tol];	
	end

	figure(1)
	hold on
	plot(1:numel(tol), (tol), 'Color', colorOrder(ibeta,:), 'LineWidth', 1.5);

	figure(2)
	hold on
	plot(1:numel(tol), log(tol), 'Color', colorOrder(ibeta,:),'LineWidth', 1.5);
	grid
end; 

figure(1)
xlabel('Iteration count');
ylabel('Error bound');
legend(legends, 'Location', 'SouthEast')
title('Error bound vs iteration count');
xlim([0 nsaiter-1+5*do_nk])

figure(2)
xlabel('Iteration count');
ylabel('log(Error bound)');
legend(legends, 'Location', 'SouthEast')
title('Error bound vs iteration count');
xlim([0 nsaiter-1+5*do_nk])
