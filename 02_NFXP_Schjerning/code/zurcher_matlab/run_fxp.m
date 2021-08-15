% run_fxp: Solve Rust's engine replacement model Rust(Ecta, 1987) 
% clear all
clc

%% Parameters for solution algorithm (used in solve.m)
ap=solve.setup;
ap.printfxp=2;						% (0= no output), (compressed output), (2= detailed iteraion output)
ap.sa_min=20;										  % Set minimum number of contraction steps (successive approximations)
ap.sa_max=ap.sa_min;	              % Set maximum number of contraction steps (successive approximations)

% algorithm switch 
algorithm = 'sa'; % sa or poly

% Read default parameters in to struct mp 
mp=zurcher.setup;


mp.integrated=0;  	% set 1 to use bellman formulated interms of 
					% integrated value function istead of expected value function	

% Transition matrix for mileage
[P] = zurcher.statetransition(mp);

% Initial guess on fixed point W
W0=0;
% bellman equation
bellman= @(V) zurcher.bellman(V, mp, P);    		
switch algorithm
	case 'sa' % solve by successive approximations (SA)
		[W, iterinfo]=solve.sa(bellman, W0, ap);	
	case 'poly' % solve using poly-algorithm (use combination of SA and NK)
		[W, pk, dV, iterinfo]=solve.poly(bellman, W0, ap, mp.beta);	
	otherwise
 		error('Algorithm must be ''sa'' or ''poly''');
end

% recall that bellman either integrated value, V or Expected value ev 
if mp.integrated;
	ev=P{1}*W;
else
	ev=W;
end


