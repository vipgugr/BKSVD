%
% BKSVD
% this function performs the Empirical Bayesian KSVD for BCD of
% histological images. Greedy version.
%
% INPUTs:
% I - RGB histological image (0-255) type double.
% K - number of columns of the given dictionary
% D0 - init D, each stain a column
%
% OUTPUTs:
% D - Color vector matrix
% X - Concentration matrix
%


function [D,X] = MB_EBKSVD4SD_v2(I,D0,K)
iter_T=0;

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
marcar=tmp>0.1; % find non-white pixels

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

term=1.e-04;
termD=5.e-03;
current_batch=0;

while current_batch<n_batches && Devol>termD
    D_batch=D;
    %Random sampling
current_batch=current_batch+1;
icol=randperm(size(Y_filtered,2),batch_size);
Y=Y_filtered(:,icol);
[P,Q] = size(Y);
X0 = D \ Y;
X0(X0 < eps) = eps;
X=X0;
minIter=2;
iter=0;
%while iter < maxIter 
while ( (iter <= minIter) || (((convH > term) || (convE > term)) && (iter <= maxIter)) )
    iter=iter+1;
      
     S_xq=(zeros(K,K,Q));

    
    parfor q = 1:Q
        pX{q} = zeros(K,1);
        pS_xq{1,1,q} = zeros(K,K);
    end
    
    parfor q=1:Q
        [xq,Sig,used,~,~,~,~,~,~] = FastLaplace(D,Y(:,q));
        pX{q}(used) = xq;
        pS_xq{1,1,q}(used,used) = Sig;
    end
    
    X = cell2mat(pX);
    %NaNs should be zeros
    X(isnan(X))=0;
    
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
        
        if K==3
            if (norm(D(:,1)-D(:,3))<0.1 || norm(D(:,2)-D(:,3))<0.1)
                D(:,3)=D0(:,3);
            end
        end
        
    end
    
    convH = sum((X(1,:)- X0(1,:)).*(X(1,:)- X0(1,:))) / sum(X0(1,:).*X0(1,:));
    convE = sum((X(2,:)- X0(2,:)).*(X(2,:)- X0(2,:))) / sum(X0(2,:).*X0(2,:));
    X0 = X;
    

    
end



Devol=norm(D_batch-D);
% D
disp(['- BKSVD - batch: ' num2str(current_batch) '- iter: ' num2str(iter) ' of ' num2str(maxIter)])

iter_T=iter_T+iter;
end

X=directDeconvolve(I,D);



end