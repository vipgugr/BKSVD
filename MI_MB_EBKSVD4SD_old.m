%
% BKSVD
% this function performs the Bayesian KSVD
%
% INPUTs:
% Images - cell object containing paths to the images
% K - number of columns of the given dictionary
% D0 - init D, each stain a column
% maxIter - maximum numer of iterations
% batch_size - num of pixel per batch
% n_batches - maximum number of batches to reach convergence
%
% OUTPUTs:
% D - dictionary



function [D] = MI_MB_BKSVD4SD(Images,D0,K,maxIter,batch_size,n_batches)
errBKSVD_T=[];
spars_T=[ ];
ne_T=[ ];

% Input check
if nargin < 5
    error('Not enough input arguments.')
end

n_images=length(Images);
images_per_batch=min(20,n_images);

i_imgs=randperm(n_images,images_per_batch);
for i=1:images_per_batch
    loaded_images{i}=double(imread(Images{i_imgs(i)}));
end
I=cell2mat(loaded_images);

Y=rgb2od(I);
[m,n,c]=size(Y);
Y_full=reshape(Y,(m)*(n),c)';

%N_side=floor(sqrt(sample_percent*m*n));
%N=N_side*N_side;
disp(['Equivalente Pixeles utilizados: ', num2str(batch_size*n_batches)])

    %Muestreo aleatorio 

icol=randperm(size(Y_full,2),batch_size);
Y=Y_full(:,icol);
%Muestreo ordenado
%Y=sort(Y_full,2,'descend');
%Y=Y(:,0.1*m*n:0.1*m*n+N);


% Initializations
[P,Q] = size(Y);
iter=0;
% Initial values
    X0 = D0 \ Y;
    X0(X0 < eps) = eps;

%D = normcols(randn(P,K));

D=D0;
D=D(:,1:K);
Devol=1;
%D =[[0.6443, 0.7167, 0.2669];[0.09, 0.9545, 0.2832];[0.6360, 0,0.7717 ]]'; 
%D= D(:,1:K);
%disp('Initializing dict with Ruifrok')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ne = zeros(1,maxIter);
spars = ne;
errBKSVD = ne;

disp(' ')
disp('Running BKSVD:')
%pbar(0,'initialize')
term=1.e-04;
current_batch=0;
while current_batch<n_batches && Devol>term
    %Muestreo aleatorio 
current_batch=current_batch+1;
% Nuevo sample of images
i_imgs=randperm(n_images,images_per_batch);
for i=1:images_per_batch
    loaded_images{i}=double(imread(Images{i_imgs(i)}));
end
I=cell2mat(loaded_images);

Y=rgb2od(I);
[m,n,c]=size(Y);
Y_full=reshape(Y,(m)*(n),c)';

tmp=mean(Y_full);
marcar=tmp>0.1; % find non-white pixels
Y_filtered=Y_full(:, marcar);

% New sample of pixels
icol=randperm(size(Y_filtered,2),batch_size);
Y=Y_filtered(:,icol);
minIter=10;
iter=0;
%while iter < maxIter 
while ( (iter <= minIter) || (((convH > term) || (convE > term)) && (iter <= maxIter)) )
    iter=iter+1;
    
    %Ultima iter para sacar la real
    %if iter==maxIter
        %Y=Y_full;
       % [P,Q] = size(Y);
    %end
    
    %initialize to zeros (otherwise at each iteration in the *used* positions new elements will be added)
%     X=(zeros(K,Q)); 
     S_xq=(zeros(K,K,Q));
% 
%     for q=1:Q
%         [xq,Sig,used,~,~,~,~,~,~] = FastLaplace_add(D,Y(:,q));
%         X(used,q)=xq;
%         S_xq(used,used,q)=Sig;
%     end
    
    parfor q = 1:Q
        pX{q} = zeros(K,1);
        %pS_xq{q} = zeros(K,K);
        pS_xq{1,1,q} = zeros(K,K);
    end
    
    parfor q=1:Q
        %[xq,Sig,used,~,~,~,~,~,~] = FastLaplace_add(D,Y(:,q));
        [xq,Sig,used,~,~,~,~,~,~] = FastLaplace(D,Y(:,q));
        pX{q}(used) = xq;
        %pS_xq{q}(used,used) = Sig;
        pS_xq{1,1,q}(used,used) = Sig;
    end
    
    X = cell2mat(pX);
    %NaNs should be zeros
    X(isnan(X))=0;
    
    S_xq = cell2mat(pS_xq);
%     parfor q = 1:Q
%         S_xq(:,:,q) = pS_xq{q};
%     end

    used_all = find(sum(abs(X),2)~=0)';
    
    if maxIter>1
        %disp('Estimando Dnew');
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
        
        Devol=norm(D-Dnew);
        %disp( norm(Devol));
        %if norm(Devol)<0.005 && iter<(maxIter-1)
        %if 
            %disp('Convergencia Dicionario')
            %iter=maxIter;
            %not_converged=false;
            
        %end
 
        D = Dnew;
        
        
    end
    
    convH = sum((X(1,:)- X0(1,:)).*(X(1,:)- X0(1,:))) / sum(X0(1,:).*X0(1,:));
    convE = sum((X(2,:)- X0(2,:)).*(X(2,:)- X0(2,:))) / sum(X0(2,:).*X0(2,:));
    X0 = X;
    
    %disp(['- BKSVD - iter: ' num2str(iter) ' of ' num2str(maxIter)])
    %pbar(iter/maxIter)
    %ne(iter)=norm(Y-D*X)/norm(Y);
    ne(iter)=Devol;
    errBKSVD(iter) = computeRMSE(D,X,Y);
    %spars(iter) = 100*(nnz(X)/numel(X));
    spars(iter)=convH+convE;
    
end
errBKSVD_T=[errBKSVD_T errBKSVD(1:iter)];
spars_T=[ spars_T spars(1:iter)];
ne_T=[ ne_T ne(1:iter)];
disp(['- BKSVD - batch: ' num2str(current_batch) '- iters: ' num2str(iter) ' of ' num2str(maxIter)])
end

X=directDeconvolve(I,D);


end