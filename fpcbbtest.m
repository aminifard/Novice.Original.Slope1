%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fpcbb %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function Output = fpcbb
%
% fbcbb is fixed point continuation method with BB steps for sparse  
%  recovery ( applied to compressed sensing and image deblurring problems)
%
%          min (1/2)*||Ax - b||^2+ mu*||x||_1 
%
% Global variables:
%
% .n        % dimension of signal vector
% .b        % a m by 1 vector
% .A        % an explicit m by n matrix, or a function handle that
%           % implements A*x and A'*x.  
% .M        % any m x m positive definite matrix, or the empty matrix.
%           % If M is empty, fpcbb assumes M = I, which reduces the first
%           % term of the objective to (1/2)*||Ax - b||_2^2 (standard
%           % least squares). 
% .options  % options for fpcbb
%         .x0            % Initial value of x
%         .mu            % regularization parameter
%         .stopCriterion % stop criterion
%         .mu0           % Initial value of mu.
%         .xs            % original signal xs. 
%         .init          % If .x0 is not defined, .init specifies how
%                        % x0 is to be initialized.  
%                           = 0 %  zeros(n,1)
%                           = 1 %  x = tau*||AtMb||_Inf * ones(n,1)
%                           = 2 %  x = tau*AtMb
%         .ftol          % tolerance on the norm(f - fp)/norm(fp)
%         .xtol          % tolerance on the norm(x - xp)/norm(xp)
%         .gtol          % tolerance on the mu*norm(g,'inf') - 1 
%         .tauD          % 0 < tauD < 2
%         .mxitr         % maximum number of the iterations 
%         .eta           % ratio of current ||b - Ax|| to approx. optimal 
%                        % ||b - Ax|| for the next mu value 
%         .fullMu        % if true, then mu = eta*sqrt(n*kap)/phi,
%                        % where phi = ||Ax-b||_M, which guarantees
%                        % that phi(now)/phi(next)>= eta. Otherwise
%                        % mu = eta*mu
%         .kappa         % required if fullMu. Is ratio of max and min 
%                        % eigenvalues of M^{1/2}AA'M^{1/2} 
%         .gamma         % weighting for reference function value, 
%                        % between 0 and 1. 0 is monotone line search,
%                        % 1 is the maximum of the N most recent values 
%                        % for the objective function 
%         .c             % constant for nonmonoton Armijo condition 
%         .beta          % step-size reduction factor for Armijo condition
%         .scale         % if true, fpcbb will calculate the maximum 
%                        % eigenvalue of A'*M*A and scale mu, A and b
%                        % appropriately if it is not equal to 1
%
% Output variables:
%
%   Output.x       % x at the last iteration
%   Output.f       % the vector of the function values for each iteration
%   Output.nf      % the number of times the function has been quantified
%   Output.itr     % the number of iterations for convergence
%                  % (or Inf if reach mxitr)
%   Output.n2re    % the vector of norm(x - xs)/norm(xs) provided that 
%                  % opts.xs exists, for each iteration.
%                  % It starts with 0th iteration.
%   Output.mse     % the value of norm(x - xs)/norm(xs) for x at last
%                  % iteration, provided that opts.xs exists.   
%   Output.cpu     % cpu time
%   Output.nz      % the number of nonzero elements in the recoverd signal 
%                  % (x at the last iteration)
%   Output.lam     % the vector of ||x||_1 for each iteration
%   Output.tau     % the vector of tau values for each iteration
%   Output.alpha   % the vector of alpha values for each iteration
%   Output.stop    % determin stop conditions
%                  % 0 % no stop
%                  % 1 % if |f-fp|<=max(1,|fp|)*ftol
%                  % 2 % if mu*||g||_inf - 1<=gtol
%                  % 3 % if ||x-xp|<=max(1,||xp||)*xtol
%                  % 4 % if i == mxiter
%
% Reference:
%
% Hale, Yin and Zhang, FPCBB.m
% 
% Written by
% Esmaeili, Shabani and Kimiaei
%


function Output = fpcbbtest

global n n_S m_S A b ub1 Mn options

mu = options.mu;
% problem dimension
m = length(b);

% unify implementation of A
if ~isa(A,'function_handle'),A = @(x,mode) explicitMatrix(A,x,mode);end


%%%%%%%%%%%%%%%%%%%%%%%%
% get or check options %
%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(options,'x0')
    if ~isvector(options.x0) || ~min(isfinite(options.x0))
        error(['If used, options.x0 should be an n x 1 vector',...
            ' of finite numbers.']);
    end
elseif isfield(options,'init')
    if ~isinInterval(options.init,0,3,true) || ...
         options.init ~= floor(options.init)
        error(['options.init must be an integer between 0 and 3,',...
            ' inclusive.']);
    end
else
    options.init = 2; 
end

if isfield(options,'mu0')
    if ~isinInterval(options.mu0,0,Inf,false)
        error('If used, options.mu0 must be a positive scalar.');
    end
end


if isfield(options,'ftol')
    if ~isinInterval(options.ftol,0,1,false)
        error(['options.ftol is tolerance on norm(x - xp)/norm(xp).',...
            ' Should be in (0,1).']);
    end
else
    options.ftol = 1E-4;
end

if isfield(options,'xtol')
    if ~isinInterval(options.xtol,0,1,false)
        error(['options.xtol is tolerance on norm(x - xp)/norm(xp).',...
           ' Should be in (0,1).']);
    end
else
    options.xtol = 1E-10;
end

if isfield(options,'gtol')
    if ~isinInterval(options.gtol,0,1,false)
        error(['options.gtol is tolerance on mu*norm(g,''inf'') - 1.',...
            ' Should be in (0,1).']);
    end
else
    options.gtol = 0.2;
end

if isfield(options,'tauD')
    if ~isinInterval(options.tauD,0,2,false)
        error('If used, options.tauD must be in (0,2).');
    end
end
    
if isfield(options,'mxitr')
    if ~isinInterval(options.mxitr,0,Inf,false) || ...
       options.mxitr~=floor(options.mxitr)
        error('options.mxitr should be a finite positive integer.');
    end
else
    options.mxitr = 10000;
end

if isfield(options,'eta')
    if ~isinInterval(options.eta,1,Inf,false)
        error('options.eta must be greater than one.');
    end
else
    options.eta = 4;
end


if isfield(options,'kappa')
    if ~isinInterval(options.kappa,1,Inf,true)
       error(['options.kappa is a condition number',...
            ' and so should be >= 1.']);
    end
end

if isfield(options,'gamma')
    if ~isinInterval(options.gamma,0,1,true)
       error('options.gamma must be in [0,1].');
    end
else
    options.gamma = 0.85;
end

if isfield(options,'c')
    if ~isinInterval(options.c,0,1,false)
       error('options.c must be in (0,1).');
    end
else
    options.c = 1E-3;
end

if isfield(options,'beta')
    if ~isinInterval(options.beta,0,1,false)
       error('options.beta must be in (0,1).');
    end
else
    options.beta = 0.5;
end

if isfield(options,'scale')
    if ~islogical(options.scale)
       error('scale should be true or false.');
    end
else
    options.scale = true;
end

%%%%%%%%%%%%%%%%%
% check scaling %
%%%%%%%%%%%%%%%%%

if options.scale
    disp(['fbcbb is checking the scaling of your',...
        ' problem because options.scale = true.']);
    eoptions.disp = 0;
    eoptions.issym = true;
    flagA = ~ isreal(A(rand(n,1),1));
    if flagA
       eoptions.isreal = false;
    end
    
   fh = @(x) A(A(x,1),2);
    
     s2     = eigs(fh,n,1,'lm',eoptions);
    flags2 = s2 > 1 - eps;
    
    if flags2
       mu = mu*s2;
       b = b/sqrt(s2);
       A = @(x,mode) A(x,mode)/sqrt(s2);
    end
    
end

% calculate A'*M*b
 AtMb = A(b,2); 

%%%%%%%%%%%%%%%%%%%%%%%%
% initialize x and tau %
%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(options,'tauD'), tauD = options.tau;
else tauD = min(1.999,-1.665*m/n + 2.665); end

if isfield(options,'x0')
    x = options.x0;
    if length(x) ~= n,error('User supplied x0 is wrong size.');end
else
    switch options.init
        case 0, x = zeros(n,1);
        case 1, x = tauD*norm(AtMb,inf)*ones(n,1);
        case 2, x = tauD*AtMb;
        case 3, x = Mn;
    end
end

%%%%%%%%%%%%%%%%%%%%%%
% stopping parameter %
%%%%%%%%%%%%%%%%%%%%%%

ftol       = options.ftol;
xtol       = options.xtol;
stopCriterion = options.stopCriterion;
mxitr      = options.mxitr;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare for iterations %
%%%%%%%%%%%%%%%%%%%%%%%%%%

Output.f     = []; Output.iterfinal  = Inf;  gp = [];  
Output.nf    = []; Output.stop= 0; Output.step = []; Output.lam = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute f,g,lam,Ax
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ax, lam
Ax = A(x,1); r = Ax - b; lam = norm(x,1);

%f
f = 0.5*(r'*r)+mu*lam;

% g
 g=A(Ax,2)-AtMb;

Output.f     = [Output.f; f];  Output.lam     = [Output.lam; lam];
nf           = 1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%
% nonmonotone and Barzilia Borwein (BB) parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%

N  = 10;lastfv(1:N)=-Inf;lastfv(N) = f; H = f; nov_k=0.85;
beta =options.beta;c=options.c;
taumin=10^-4;taumax=10^4;warned =false;

tic;
for i = 1:mxitr
    
    % get bb step
    flagnf = nf > 1;
    if flagnf
        dg = g - gp;
        tau = nrmxxp/max(real(xxp'*dg),eps);
        tau = min(taumax,max(taumin,tau));
    else tau = tauD;
    end
    nu = tau*mu  ;
    
    % store old point
    xp = x; fp=f; gp = g; Axp = Ax;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % take fixed-point step and compute shrinkage direction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y = x - tau*g;
    x = sign(y).*max(0,abs(y)-nu);
    x(n_S+1:n)=max(0,x(n_S+1:n));
    x(1:n_S)=ub1-x(n_S+1:2*n_S);
    dx = x - xp;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % calculate g, f, Ax, lam
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Ax, lam
    Ax = A(x,1);r = Ax - b; lam = norm(x,1);
    %f
    f =0.5*(r'*r)+mu*lam;
    nf = nf+1;
    
    %  g
    g=A(Ax,2)-AtMb;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Nonmonotone line search %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
      dg = g - gp; dAx = Ax - Axp;
      
    const = c*((mu*sign(xp) + gp)'*dx);
    if ~isreal(const) && ~warned
        warning('FPC-BB is designed for problems with real data.'); 
        warned = true;
    end
    % xnew = xp + alpha*dx
     alpha = 1;  cnt = 0;
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while 1
         flagarmijo = f <= H + alpha*const;
         flagcnt = cnt==100;
         if flagarmijo|flagcnt, break;end
         % update alpha, x, Ax, lam,f and g
         alpha = alpha*beta; x = xp + alpha*dx; 
         x(n_S+1:n)=max(0,x(n_S+1:n));
         x(1:n_S)=ub1-x(n_S+1:2*n_S);
         Ax    = Axp + alpha*dAx; r = Ax - b;lam = norm(x,1);
         f = 0.5*norm(r)^2 + mu*lam;
         nf = nf+1; g = gp + alpha*dg; 
         cnt = cnt +1;
    end
      
    %%%%%%%%%%%%%%%%%%%%%%%%
    % Nonmonotone Strategy %
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    modiN  = mod(i,N);
    flagNM = (modiN== 0);
    if flagNM, lastfv(N) = f;else lastfv(modiN) = f;end
    
    fmax = max(lastfv); Ng   = norm(g);
    
    
    flagNg = Ng<=0.01;
    if flagNg,nov_k=(2/3)*nov_k+0.01;else nov_k= max(0.99*nov_k,0.5);end
   
    H    = nov_k*fmax+(1-nov_k)*f;
   %%%%%%%%%%%%%%
    % Outputs %
    %%%%%%%%%%%%%% 
    Output.f = [Output.f; f];  Output.lam  = [Output.lam; lam]; 
    ffp = f - fp; nrmffp   = norm(ffp); xxp  = x - xp; nrmxxp  = norm(xxp);
    Output.step  = [Output.step; nrmxxp]; 
    
    %%%%%%%%%%%%%%%%%%
    % Stop condition %
    %%%%%%%%%%%%%%%%%%
    
    switch stopCriterion
        
        case 1
            
            crit1     = nrmffp/max(norm(fp),1); 
            flagcrit1 = (crit1 < ftol);
            if flagcrit1, Output.stop=1;end
            
       
        case 2
            
            crit2     = nrmxxp/max(norm(xp),1); 
            flagcrit2 = (crit2 < xtol);
            if flagcrit2, Output.stop=3;end
            
        case 3
            
           if i==mxitr, Output.stop=4;end
    end
     
    x(n_S+1:n)=max(0,x(n_S+1:n))
    x(1:n_S)=ub1-x(n_S+1:2*n_S)
    if ismember(Output.stop,[1 2 3]) | i==mxitr
       e = toc; nz_x = (x ~= 0.0);num_nz_x = sum(nz_x(:));
       Output.cpu = e; Output.x= x; Output.iterfinal = i; Output.nf  = nf; 
       Output.nz = num_nz_x ; 
       break;
    end

end
end % fpcbb

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%