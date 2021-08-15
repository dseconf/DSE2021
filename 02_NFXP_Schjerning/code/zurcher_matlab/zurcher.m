classdef zurcher
    % zurcher class: Holds model parts for Rust's engine repplacement model Rust(Ecta, 1987) 
    % By Fedor Iskhakov, Bertel Schjerning, and John Rust
    methods (Static)
        function mp = setup(mpopt)
            % zurcher.setup: Sets default params equal to parameters in Table X in Rust (1987).    
            % Spaces
            mp.n=175;                       % Number of grid points     
            mp.max=450;                     % Max of mileage

            % Structural parameters
            mp.p=[0.0937 0.4475 0.4459 0.0127]';    % Transition probabilities
            mp.RC=11.7257;                          % Replacement cost
            mp.c=2.45569;                           % Cost parameter
            mp.beta=0.9999;                         % Discount factor

            mp.integrated=0;                        % use integrated value function instead 
                                                    % of expected value function
            if nargin>0
              pfields=fieldnames(mpopt);
              for i=1:numel(pfields);
                  mp.(pfields{i})=mpopt.(pfields{i});
              end
            end
        
            mp.grid= (0:mp.n-1)';           % Grid over mileage
        end % end of zurcher.setup

        function P = statetransition(mp)
            % zurcher.statetransition: Computes state transition matrices conditional on choice. 
            %
            % Inputs
            %  mp       structure with model parameters (see setup)
            %
            % Outputs:
            %  P        2 dimensional cell array with mp.n x mp.n conditional transition matrices
            %  where    P{1} is transition matrix conditional on keeping engine
            %  and      P{2} IS Transition matrix conditional on replacing engine

            % Transition matrix conditional on keeping engine  
            p=[mp.p; (1-sum(mp.p))];
            n=mp.n;
            P{1}=0;
            for i=0:numel(p)-1;
                P{1}=P{1}+sparse(1:n-i,1+i:n,ones(1,n-i)*p(i+1), n,n);
                P{1}(n-i,n)=1-sum(p(1:i));
            end
            P{1}=sparse(P{1});

            % transition conditional on replacing engine
            P{2}=sparse(n,n); 
            for i=1:numel(p);
                P{2}(:,i)=p(i);
            end
        end % end of zurcher.statetransition

        function [W1, pk, dBellman_dV]=bellman(W0, mp, P)
            % zurcher.bellman:   Procedure to compute Bellman operator
            % If mp.integrated == 1  W1 and W0 is the integrated value function, V1=W1 given V0=W0. 
            % If mp.integrated == 0  W1 and W0 is the expected value function, ev1=W1 given ev0=W0
            % To obtain ev=bellamn_ev(ev) from V=bellamn_integrated(V) compute ev=P{1}*V

            % Inputs and outputs: see bellman_integrated and bellman_ev
            if mp.integrated
                [W1, pk, dBellman_dV]=zurcher.bellman_integrated(W0, mp, P);
            else
                [W1, pk, dBellman_dV]=zurcher.bellman_ev(W0, mp, P);
            end
        end % end of zurcher.bellman

        function [V1, pk, dBellman_dV]=bellman_integrated(V0, mp, P)
            % zurcher.bellman_integrated:    Procedure to compute bellman integrated Bellman equation
            %
            % Inputs
            %  V0       mp.n x 1 matrix of initial guess on integrated value function
            %  mp       structure with model parameters (see setup)
            %  P        2 dimensional cell array with mp.n x mp.n conditional transition matrices (see statetransition)
            %
            % Outputs:
            %  V1       mp.n x 1 matrix of integrated value function given initial guess V0
            %  pk       mp.n x 1 matrix of choice probabilities (Probability of keep)

            if nargin<3
                P = zurcher.statetransition(mp);
            end

            if numel(V0)==1;
                V0=zeros(mp.n,1);
            end

            cost=0.001*mp.c*mp.grid;                    % Cost function
            vK=-cost            + mp.beta*P{1}*V0;      % Value off keep
            vR=-mp.RC-cost(1)   + mp.beta*P{2}*V0;      % Value of replacing
            maxV=max(vK, vR);
            V1=(maxV + log(exp(vK-maxV)  +  exp(vR-maxV))); 

            % If requested, also compute choice probability from ev (initial input)
            if nargout>1 
                pk=1./(1+exp((vR-vK)));
            end

            if nargout>2 % compute Frechet derivative
                dBellman_dV=mp.beta*(bsxfun(@times,P{1},pk) + bsxfun(@times,P{2},1-pk));
            end
        end % end of zurcher.bellman_integrated  

        function [ev1, pk, dbellman_dev]=bellman_ev(ev, mp, P)
            % zurcher.bellman_EV: Procedure to compute bellman equation based on expected values
            %                     (Similar to Rust NFXP_Manual)
            % Inputs
            %  ev       mp.n x 1 matrix of expected values given initial guess on value function
            %  mp       structure with model parameters (see setup)
            %  P        2 dimensional cell array with mp.n x mp.n conditional transition matrices (see statetransition)
            %
            % Outputs:
            %  ev1      mp.n x 1 matrix of expected values given initial guess of ev 
            %  pk       mp.n x 1 matrix of choice probabilities (Probability of keep)

            cost=0.001*mp.c*mp.grid;             % Cost function

            vK=-cost          + mp.beta*ev;      % Value off keep
            vR=-cost(1)-mp.RC + mp.beta*ev(1);   % Value of replacing    

            % Need to recentered logsum by subtracting max(vK, vR)
            maxV=max(vK, vR);
            V=(maxV + log(exp(vK-maxV)  +  exp(vR-maxV))); 
            ev1=P{1}*V;

            % If requested, also compute choice probability from ev (initial input)
            if nargout>1 
                pk=1./(1+exp((vR-vK)));
            end
            if nargout>2 % compute Frechet derivative
                dbellman_dev=sparse(mp.n,mp.n);
                dbellman_dev=mp.beta*bsxfun(@times, P{1}, pk');    
                dbellman_dev(:,1)=dbellman_dev(:,1)+mp.beta*P{1}*(1-pk);        % Add additional term for derivative wrt Ev(1), since Ev(1) enter logsum for all states            
            end
        end % end of zurcher.bellman_ev  

        function [f,g,h]=ll_integrated(data, mp, pnames, theta, ap)
            global V0;

            if ~exist('ap') 
                ap=solve.setup;
            end

            % update model parameters
            mp=vec2struct(theta, pnames, mp);

            mp.p=abs(mp.p); %helps BHHH which is run as unconstrained optimization
            n_c=numel(mp.c);
            N=size(data.x,1);

            % Update P 
            P = zurcher.statetransition(mp);            

            % Solve model
            bellman= @(V) zurcher.bellman_integrated(V, mp, P);
            [V0, pk, dV, iterinfo]=solve.poly(bellman, V0, ap, mp.beta);  
            F=speye(mp.n) - dV;
            
            % Evaluate likelihood function
            % log likelihood regarding replacement choice
            pxK=pk(data.x);
            pxR=1-pxK;

            logl=log(pxK.*(1-data.d) + pxR.*data.d);

            % add on log like for mileage process
            if numel(theta)>2
                p=[mp.p; 1-sum(mp.p)];
                n_p=numel(p)-1;
                logl=logl + log(p(1+ data.dx1));
            else
                n_p=0;
            end

            % Objective function (negative mean log likleihood)
            f=mean(-logl);

            if nargout >=2;  %% compute scores
                dc=0.001*mp.grid;
                cost=0.001*mp.c*mp.grid;

                % step 1: compute derivative of contraction operator wrt. parameters
                dbellman_dmp=   zeros(mp.n,1+ n_c + n_p);
                dbellman_dmp(:, 1)= P{1}*pk-1;                  % Derivative wrt. RC
                dbellman_dmp(:, 2:1+n_c)=-(P{1}*dc).*pk;       % Derivative wrt. c

                if numel(theta)>2
                    vk=-cost+mp.beta*P{1}*V0;
                    vr=-mp.RC-cost(1)+mp.beta*P{2}*V0;
                    vmax=max(vk,vr);
                    dbellman_dPi=vmax+log(exp(vk-vmax)+exp(vr-vmax));

                    for iP=1:n_p;
                        dbellman_dmp(1:mp.n-iP, 1+n_c+iP)= [dbellman_dPi(iP:end-1)] - [dbellman_dPi(n_p+1:mp.n); repmat(dbellman_dPi(end), n_p-iP, 1)];
                    end
                    invp=exp(-log(p));
                    invp=[sparse(1:n_p,1:n_p,invp(1:n_p),n_p,n_p); -ones(1,n_p)*invp(n_p+1)];
                    N=size(data.x,1);
                end
                
                % step 2: compute derivative of ev wrt. parameters
                dVdmp=F\dbellman_dmp;  
                
                % step 3: compute derivative of log-likelihood wrt. parameters
                score=bsxfun(@times, (data.d-pxR),[-ones(N,1) dc(data.x,:) zeros(N,n_p)] + (dVdmp(ones(N,1),:)-dVdmp(data.x,:)));    
                if numel(theta)>2
                    for iP=1:n_p;
                        score(:,1+n_c+iP)= score(:,1+n_c+iP) + invp(1+data.dx1,iP);
                    end
                end
                g=mean(-score,1);
            end

            if nargout >=3;  %% compute Hessian
                h=score'*score/(size(logl, 1)); 
            end 
        end % end of zurcher.ll_integrated 
                                
        function [f,g,h]=ll_ev(data, mp, pnames, theta, ap)

            global ev0;

            if ~exist('ap') 
                ap=solve.setup;
            end


            % update model parameters
            mp=vec2struct(theta, pnames, mp);

            mp.p=abs(mp.p); %helps BHHH which is run as unconstrained optimization
            n_c=numel(mp.c);
            N=size(data.x,1);

            % Update P 
            P = zurcher.statetransition(mp);


            % Solve model
            bellman= @(ev) zurcher.bellman_ev(ev, mp, P);
            [ev0, pk, dev, iterinfo]=solve.poly(bellman, ev0, ap, mp.beta);  
            F=speye(mp.n) - dev;
            
            % Evaluate likelihood function
            % log likelihood regarding replacement choice
            pxK=pk(data.x);
            pxR=1-pxK;

            logl=log(pxK.*(1-data.d) + pxR.*data.d);

            % add on log like for mileage process
            if numel(theta)>2
                p=[mp.p; 1-sum(mp.p)];
                n_p=numel(p)-1;
                logl=logl + log(p(1+ data.dx1));
            else
                n_p=0;
            end

            % Objective function (negative mean log likleihood)
            f=mean(-logl);

            if nargout >=2;  %% compute scores
                dc=0.001*mp.grid;
                cost=0.001*mp.c*mp.grid;

                % step 1: compute derivative of contraction operator wrt. parameters
                dbellman_dmp=zeros(mp.n,1+ n_c + n_p);

                dbellman_dmp(:, 1)=(1-pk)*(-1);         % Derivative wrt. RC
                dbellman_dmp(:, 2:1+n_c)=pk.*(-dc);     % Derivative wrt. c

                if numel(theta)>2
                    vk=-cost+mp.beta*ev0;
                    vr=-mp.RC-cost(1)+mp.beta*ev0(1);
                    vmax=max(vk,vr);
                    dbellman_dPi=vmax+log(exp(vk-vmax)+exp(vr-vmax));

                    for iP=1:n_p;
                        dbellman_dmp(1:mp.n-iP, 1+n_c+iP)= [dbellman_dPi(iP:end-1)] - [dbellman_dPi(n_p+1:mp.n); repmat(dbellman_dPi(end), n_p-iP, 1)];
                    end
                    invp=exp(-log(p));
                    invp=[sparse(1:n_p,1:n_p,invp(1:n_p),n_p,n_p); -ones(1,n_p)*invp(n_p+1)];
                    N=size(data.x,1);
                end
                
                % step 2: compute derivative of ev wrt. parameters
                devdmp=F\dbellman_dmp;  
                
                % step 3: compute derivative of log-likelihood wrt. parameters
                score=bsxfun(@times, (data.d-pxR),[-ones(N,1) dc(data.x,:) zeros(N,n_p)] + (devdmp(ones(N,1),:)-devdmp(data.x,:)));    
                if numel(theta)>2
                    for iP=1:n_p;
                        score(:,1+n_c+iP)= score(:,1+n_c+iP) + invp(1+data.dx1,iP);
                    end
                end
                g=mean(-score,1);
            end

            if nargout >=3;  %% compute Hessian
                h=score'*score/(size(logl, 1)); 
            end 
        end % end of zurcher.ll_ev

        function [data] = simdata(N,T,mp, P, pk, seed)
            % zurcher.simdata: simulates data from engine replacement model. 
            %
            % Inputs
            %  N       Number busses to simulate
            %  T       Number of time periods to be simulated for each bus
            %  mp       structure with model parameters (see setup)
            %  P        2 dimensional cell array with mp.n x mp.n conditional transition matrices (see statetransition)
            %  pk       mp.n x 1 matrix of choice probabilities (probability of keep)
 
            % Initial conditions: all buses start uniformly distributed across statespace

            % Outputs:
            %  data: table with NxT rows and 6 columns 
            %       data.id  :  Bus id 
            %       data.t   :  Time period
            %       data.d   :  Replacement indicator (Replace = 1, Keep=0)
            %       data.x   :  Mileage in period t
            %       data.x1  :  Mileage in period t+1
            %       data.dx1 :  Change in mileage between period t and t+1

            if exist('seed')
                rng(seed,'twister');
            end

            id=repmat((1:N),T,1);
            t=repmat((1:T)',1,N);
            u_init=randi(mp.n,1,N);
            u_dx=rand(T,N);
            u_d=rand(T,N);
            
            csum_p=cumsum(mp.p);
            dx1=0;
            for i=1:numel(csum_p);
                dx1=dx1+ (u_dx>(csum_p(i)));
            end;
            
            x=nan(T, N);
            x1=nan(T, N);
            x(1,:)=u_init; % Intial conditions
            
            for it=1:T;
                d(it,:)=(u_d(it,:)<(1-pk(x(it,:)')'))*1;  % Replace = 1, Keep=0
                % x1(it,:)=min(x(it,:).*(d(it,:)==1) +(d(it,:)==2) + dx1(it,:), mp.n);
                x1(it,:)=min(x(it,:).*(1-d(it,:)) + d(it,:) + dx1(it,:), mp.n);
                if it<T;
                    x(it+1,:)=x1(it,:);
                end
            end
            
            data.id=id;
            data.t=t;
            data.d=d;
            data.x=x;
            data.x1=x1;
            data.dx1=dx1;
            
            pfields=fieldnames(data);
            for i=1:numel(pfields);
                data.(pfields{i})=reshape(data.(pfields{i}), T*N,1);
            end
            data=struct2table(data);
        end % end of zurcher.simdata
        
        function dta = readbusdata(mp,bustypes)
            load('busdata1234.mat');

            % Select busses
            if nargin>1
                data=data(data(:,2) <=bustypes,:);
            end

            id=data(:,1);       % Bus id
            bustype=data(:,2);  % bustype: 1,2,3,4
            d1=data(:,5);       % Lagged replacement dummy, d_(t-1)
            d=[d1(2:end);0];    % Replacement dummy, d_t
            x=data(:,7);        % Odometer, x_t


            % Discretize odometer data into 1, 2, ..., n
            x=ceil(x.*mp.n/(mp.max*1000));                                 

            % Mothly milage
            dx1=x(1:end,1)-[0;x(1:end-1,1)];           % Take first difference on odometer                           
            dx1=dx1.*(1-d1)+x.*d1;                     % Make true x_t-x_(t-1) (replace first difference by x_t if replacement dummy lagged is 1)

            % remove observations with missing lagged mileage
            data=[id, bustype, d,x,dx1];                                           % [ id , replace dummy lagged, odometer, "change in odometer" , bustype]
            remove_first_row_index=data(:,1)-[0;data(1:end-1,1)];                  % Take first difference of ID...
            data=data((remove_first_row_index==0),:);                              % ... and only keep lines where ID hasn't changed 

            % Save data structure
            dta.d=data(:,3);           
            dta.x=data(:,4);
            dta.dx1=data(:,5);
        end % end of zurcher.readbusdata    

        function [pp, pp_K, pp_R] = eqb(mp,P, pk)
            % This program computes the equilibrium distribution of the controlled
            % stochastic process P(i(t)|x(t))p(x(t)|x(t-1),i(t-1)). The equilibrium
            % distribution pp is an (1,n) matrix, where n is the number of discrete
            % mileage cells: pp=pp*P, where P is the transition probability matrix
            % for the (nxn) merged state process, where states with i=1 are 
            % identified with the state x=1. Estimates of the complete (1,2*n)
            % equilibrium distribution are then uncovered from pp. 

            % Outputs    
            % pp: Pr{x} (Equilibrium distribution of mileage)
            % pp_K: Pr{x,i=Keep}
            % pp_R: Pr{x,i=Replace}

            function ed=ergodic(p);
                % ergodic.m: finds the invariant distribution for 
                % an NxN Markov transition probability
                n=size(p,1);
                if n ~= size(p,1);
                    fprint('Error: p must be a square matrix\n');
                    ed=NaN;
                else
                    ap=eye(n)-p';
                    ap=[[ap; ones(1,n)], ones(n+1,1)];
                    if (rank(ap) < n+1);
                        fprintf('Error: transition matrix p is not ergodic\n');
                        ed=NaN;
                    else
                        ed=[ones(n,1); 2];
                        ed=ap\ed;
                        ed(end)=[];
                        ed=reshape(ed,1,[]); %row vector
                    end;
                end
            end % end of ergodic

            % state transition matrix 
            pl = P{1}.*pk +P{2}.*(1-pk);

            % Equilibrium distribution of mileage
            pp = ergodic(pl);   % Pr{x} (Equilibrium distribution of mileage)
            pp_K=pp.*(pk');     % Pr{x,i=Keep}
            pp_R=pp.*(1-pk');   % Pr{x,i=Replace}

        end % end of zurcher.eqb
  
    end % end of methods
end % end of zurcher class