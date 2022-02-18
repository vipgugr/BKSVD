%
% BKSVD
% this function performs the Bayesian KSVD for Blind Color Deconvolution of
% Histological images. Full Inference.
%
% INPUTs:
% I - RGB histological image (0-255) type double.
% K - number of stains (usually K=2)
% D0 - init D, each stain a column (Ruifrok's reference stain matrix)
%
% OUTPUTs:
% D - Color Vector Matrix
% X - Concentration matrix
% 
%
%


function [D,X] = MB_BKSVD4SD_FI_v3(I,D0,K)
iter_T=0;
% Input check
if nargin < 3
   error('Not enough input arguments.')
end
    batch_size=1000;
    n_batches=10;
    maxIter=100;

Y=rgb2od(I);
[m,n,c]=size(Y);
Y_full=reshape(Y,(m)*(n),c)';

tmp=mean(Y_full);
marcar=tmp>0.2; % find non-white pixels

Y_filtered=Y_full(:, marcar);
if size(Y_filtered,2)<batch_size
	disp('batch_size reduced, not enough pixels')
	batch_size=size(Y_filtered,2)
    n_batches=1
end

D=D0;
D=D(:,1:K);
Devol=1;
%D =[[0.6443, 0.7167, 0.2669];[0.09, 0.9545, 0.2832];[0.6360, 0,0.7717 ]]'; 
%D= D(:,1:K);
%disp('Initializing dict with Ruifrok')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp(' ')
disp('Running BKSVD:')


%Init sigma noise
sigma_noise =  std(Y_full(:))^2/1e2;
a_noise=1;
b_noise=0.5;
term=1.e-04;
termD=5.e-03;
current_batch=0;
while current_batch<n_batches && Devol>termD
    %Random sampling
    D_batch=D;
    current_batch=current_batch+1;
    icol=randperm(size(Y_filtered,2),batch_size);
    Y=Y_filtered(:,icol);
    [P,Q] = size(Y);
    X0 = D \ Y;
    X0(X0 < eps) = eps;
    X=X0;

    %inits

    gamma=(ones(K,Q));
    S_xq=(zeros(K,K,Q));
    minIter=2;
    iter=0;
    %while iter < maxIter 
    while ( (iter <= minIter) || (((convH > term) || (convE > term)) && (iter <= maxIter)) )
    iter=iter+1;
%     
    parfor q = 1:Q
        pX{q} = zeros(K,1);
        p_gamma{q} = zeros(K,1);
        %pS_xq{q} = zeros(K,K);
        pS_xq{1,1,q} = zeros(K,K);
    end
    
    parfor q=1:Q
        warning('off')
		
		%Calcular lambda
		lambda_q= (2*K-2) / sum(gamma(:,q));
        %Calcular gamma
        
        for s=1:K
            gamma_s1=-1/(2*lambda_q) + sqrt(1/(4*lambda_q^2) + (X0(s,q)^2+S_xq(s,s,q))/lambda_q);
            if gamma_s1==0 gamma_s1=eps; end
            p_gamma{q}(s)=gamma_s1;
        end
        
        %Calcular S_x
        pS_xq{1,1,q}=(sigma_noise*D'*D + diag(p_gamma{q})^-1)^-1;
        %Calcular x
        pX{q} =sigma_noise*pS_xq{1,1,q}*D'*Y(:,q);     
        w = warning('query','last')  

		
    end
    

    X = cell2mat(pX);
    %NaNs should be zeros
    X(isnan(X))=0;
    gamma = cell2mat(p_gamma);
    S_xq = cell2mat(pS_xq);

    used_all = find(sum(abs(X),2)~=0)';
    
    if maxIter>1

        Dnew = D; % Not to mix new and old elements in D while updating the dictionary

        % estimation of D
        ak=zeros(P,K); bk=ak; ck=zeros(1,K); ek=ck; tk=ak;
        Sq=sum(S_xq,3);
        for z = 1:numel(used_all)
            k = used_all(z);
            ak(:,k)=sum(D(:,[1:(k-1) (k+1):K])*Sq([1:(k-1) (k+1):K],k),2);
            bk(:,k)=(Y-D*X+D(:,k)*X(k,:))*X(k,:)';
            ck(k)= sum(S_xq(k,k,:));
            ek(k)=sum(X(k,:).^2)+ck(k);
            tk(:,k)=1/sqrt(ek(k))*(bk(:,k)-ak(:,k));
            Dnew(:,k)=tk(:,k)/norm(tk(:,k));
        end
         
        D = Dnew;
        
        %This avoid the third vector to replace the others
        if K==3
            if (norm(D(:,1)-D(:,3))<0.1 || norm(D(:,2)-D(:,3))<0.1)
                D(:,3)=D0(:,3);
            end
        end
        
    end
    
    %Restimate noise
    sigma_noise= (3*Q+2*a_noise)/ (norm(Y-D*X)^2+2*b_noise);
    
    convH = sum((X(1,:)- X0(1,:)).*(X(1,:)- X0(1,:))) / sum(X0(1,:).*X0(1,:));
    convE = sum((X(2,:)- X0(2,:)).*(X(2,:)- X0(2,:))) / sum(X0(2,:).*X0(2,:));
    X0 = X;
    

    
end
Devol=norm(D_batch-D);
disp(['- BKSVD - batch: ' num2str(current_batch) '- iter: ' num2str(iter) ' of ' num2str(maxIter)])
iter_T=iter_T+iter;% Devol_T=[Devol_T Devol];
end

X=directDeconvolve(I,D);
%X(X < eps) = eps;


end