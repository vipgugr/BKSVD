function [weights,Sig,used,sigma2,errbars,basis,selected,alpha,lambda] = FastLaplaceAddSigmaIn(PHI,y,sigma2,eta,lambda_init, adaptive,optimal,scale, verbose)
% This code implements the fast Laplace algorithm from the following paper:
% [1] S. D. Babacan, R. Molina, A. K. Katsaggelos. �Bayesian Compressive Sensing using Laplace Priors,�
% submitted for publication, IEEE Transactions on Image Processing, September 2008.
%
% This code is based on the BCS code available from http://people.ee.duke.edu/~lihan/cs/ from the following paper
% [2] S. Ji, Y. Xue, L. Carin, "Bayesian Compressive Sensing," IEEE Trans. Signal Processing, vol. 56, no. 6, June 2008.
%
% Please check the accompanying license and the license of [2] before using.
% Note that this code is not optimized for speed, and it is generalized to investigate different choices of the
% parameters.
%
% Author (aka person to blame): Derin Babacan, sdb@northwestern.edu
% Last modified by: 08/08/08
%
% Usage:
% [weights,used,sigma2,errbars,basis,alpha,lambda] = FastLaplace(PHI,y,sigma2,eta, lambda_init, adaptive,optimal,scale, verbose)
% Inputs:
%   Required:
%       PHI: measurement matrix
%       y:   CS measurements
%   Optional:
%       sigma2: initial noise variance (default : std(t)^2/1e2)
%       eta:  threshold for stopping the algorithm (default : 1e-8)
%       lambda_init : To set lambda equal to a fixed nonegative value.
%                     if lambda_init = [], lambda will be computed automatically, which is the suggested method.
%                     lambda_init = 0 corresponds to the BCS algorithm in [2], see [1] for technical details.
%
%   Inputs for Adaptive CS (this part is left unchanged from the BCS code, see [2])
%       adaptive: generate basis for adpative CS? (default: 0)
%       optimal: use the rigorous implementation of adaptive CS? (default: 1)
%       scale: diagonal loading parameter (default: 0.1)
% Outputs:
%   weights:  sparse weights
%   used:     the positions of sparse weights
%   sigma2:   re-estimated noise variance
%   errbars:  one standard deviation around the sparse weights
%   basis:    if adaptive==1, then basis = the next projection vector, see [2]
%   alpha:    sparse hyperparameters (1/gamma), see [1]
%   lambda:   parameter controlling the sparsity , see [1]

%% Check inputs
if nargin < 2,
    error('Not enough inputs');
end
if nargin < 3,
    sigma2 =  std(y)^2/1e2;
end
if nargin < 4,
    eta = 1e-8;
end
if nargin < 5,
    lambda_init = [];
end
if nargin < 6
    adaptive = 0;
end
if nargin < 7
    optimal = 1;
end
if nargin < 8
    scale = 0.1;
end
if nargin < 9,
    verbose = 0;
end

% find initial alpha
[N,M] = size(PHI);
PHIy = PHI'*y;
PHI2 = sum(PHI.^2)';
ratio = (PHIy.^2)./PHI2;
[maxr,index] = max(ratio);

alpha = PHI2(index)/(maxr-sigma2);
% compute initial mu, Sig, S, Q
phi = PHI(:,index);
Hessian = alpha + phi'*phi/sigma2;
Sig = 1/Hessian;
mu = Sig*PHIy(index)/sigma2;
left = PHI'*phi/sigma2;
S = PHI2/sigma2-Sig*left.^2;
Q = PHIy/sigma2-Sig*PHIy(index)/sigma2*left;

% Keep track of the positions selected during the iterations
selected = index;
max_it = 100;

nu = 0;

membership = zeros(1,M);
membership(index) = 1;

for count = 1:max_it
    
    s = S; q = Q; %this assumes that all gamma_i = 0
    s(index) = alpha.*S(index)./(alpha-S(index)); % and just changes one (some) where alpha_i~=0
    q(index) = alpha.*Q(index)./(alpha-S(index));
    
    if isempty(lambda_init),
        lambda = (nu + 2*( length(index) - 1 )) / (nu + sum(1./alpha));
    else
        lambda = lambda_init;
    end
    
    
    A = lambda + s - q.^2;
    B = 2*lambda.*s + s.^2;
    C = lambda.*s.^2;
    
    theta = q.^2-s;
    
    discriminant = B.^2 - 4.*A.*C;
    
    nextAlphas = (-B - sqrt(discriminant) ) ./ (2*A);
    
    % choose the next alpha that maximizes marginal likelihood
    ml = -inf*ones(1,M);
    
    update = false;
    for m=1:M
        if ~membership(m)
            if theta(m) > lambda
                Alpha = nextAlphas(m);
                ml(m) = log(Alpha ./ (Alpha + s(m)) )+ q(m).^2 ./ (Alpha + s(m)) - lambda./Alpha;
                update = true;
            end
        end
    end
    
    if ~update
        break
    end
    
    [ML(count),idx] = max(ml);
    
    % check convergence
    if count > 2 && abs(ML(count)-ML(count-1)) < abs(ML(count)-ML(1))*eta
        break;
    end
    
    % update alphas
    % Choose the basis which results in the largest increase in the
    % likelihood

    
    if theta(idx) > lambda
        if isempty(find(index==idx,1)) % reestimate a basis
            Alpha = nextAlphas(idx);
            phii = PHI(:,idx); Sigii = 1/(Alpha+S(idx));
            mui = Sigii*Q(idx);
            comm1 = Sig*(phi'*phii)/sigma2;
            ei = phii-phi*comm1;
            off = -Sigii*comm1;
            %
            Sig = [Sig+Sigii*(comm1*comm1'), off; off', Sigii];
            mu = [mu-mui*comm1; mui];
            comm2 = PHI'*ei/sigma2;
            S = S - Sigii*comm2.^2;
            Q = Q - mui*comm2;
            %
            index = [index;idx]; membership(idx) = 1;
            %
            alpha = [alpha;Alpha];
            phi = [phi,phii];           
        end
    end
    
    selected = [selected; idx];
    
end

weights	= mu;
used = index;
errbars = sqrt(diag(Sig));

basis = [];
% generate a basis for adaptive CS?
if adaptive
    if optimal
        [V,D] = eig(Sig);
        [~,idx] = max(diag(D));
        basis = V(:,idx)';
    else
        temp = phi'*phi/sigma2;
        Sig_inv = temp + scale*mean(diag(temp))*eye(length(used));
        [V,D] = eig(Sig_inv);
        [~,idx] = min(diag(D));
        basis = V(:,idx)';
    end
end
%fprintf(1,'Algorithm converged, # iterations : %d \n',count);

