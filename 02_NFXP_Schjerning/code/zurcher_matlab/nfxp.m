classdef nfxp
  % NFXP class: Estimates rust's engine repplacement model Rust(Ecta, 1987) 
  % By Fedor Iskhakov, Bertel Schjerning, and John Rust
  methods (Static)
    function [results, pnames, theta_hat, Avar] = estim(data, mp, twostep)
      samplesize=numel(data.d);
 
      % set options optimizer options for fminunc      
      optim_options= optimset('Algorithm','trust-region', 'Display','off', 'GradObj','on', 'TolFun',1E-7,'TolX',1E-7,'Hessian','on');

      % ************************************
      % STEP 1: ESTIMATE p 
      % ************************************

      % Estimate model using Partial MLE
      tab=tabulate(data.dx1);
      tab=tab(tab(:,3)>0,:);
      p=tab(1:end-1,3)/100;

      % ************************************
      % STEP 2: ESTIMATE structural parameters
      % ************************************

      % Use first step estimates as starting values for p
      mp.p=p;

      outsidetimer=tic;
      
      % Estimate RC and c
      pnames={'RC', 'c'};
      if mp.integrated
        llfun=@(theta) zurcher.ll_integrated(data, mp, pnames, theta);  
      else
        llfun=@(theta) zurcher.ll_ev(data, mp, pnames, theta);  
      end

      theta_start=struct2vec(mp,pnames);
      [theta_hat,FVAL,EXITFLAG,OUTPUT1]=fminunc(llfun ,theta_start,optim_options);

      % Estimate RC, c and p
      if twostep==0;
        pnames={'RC', 'c', 'p'};
        if mp.integrated
          llfun=@(theta) zurcher.ll_integrated(data, mp, pnames, theta);  
        else
          llfun=@(theta) zurcher.ll_ev(data, mp, pnames, theta); 
        end
        theta_start=struct2vec(mp,pnames);
        [theta_hat,FVAL,EXITFLAG,OUTPUT2]=fminunc(llfun,theta_start,optim_options);
      end

      timetoestimate=toc(outsidetimer);

      % Compute Variance-Covaiance matrix
      [f,g,h]=llfun(theta_hat);
      Avar=inv(h*samplesize);
      
      % Collect results
      results=vec2struct(theta_hat,{'RC', 'c'},mp);
      results.cputime=timetoestimate;
      results.converged   =   (EXITFLAG>=1 &  EXITFLAG<=3);

      results.MajorIter=OUTPUT1.iterations;
      results.funcCount=OUTPUT1.funcCount;
      if twostep==0
        results.MajorIter=OUTPUT1.iterations+OUTPUT2.iterations;
        results.funcCount=OUTPUT1.funcCount+OUTPUT2.funcCount;
      end
      results.llval   =   -FVAL*samplesize;
    end % end of nfxp.estim
  end % end of methods
    
    
end % end of estim class