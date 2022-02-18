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



function [D, C_Rmax] = MI_MB_BKSVD4SD(Images,D0,K)
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

[D,CT] = MB_EBKSVD4SD(double(I),D0,K);
C_Rmax = prctile(CT',99)


end