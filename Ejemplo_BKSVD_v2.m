%% Params
clc,clear all, close all
ns=2;
maxIter=100;
batch_size=1000;
n_batches=10;
%Initial reference matrix
load 'MLandini' RM;
D0=RM(:,1:ns);

%% Reference Image


I_ref=imread('Reference.jpg');
[Mref,Cref] = MB_BKSVD4SD(double(I_ref),D0,ns);
Cref_Rmax = prctile(Cref',99);
%%
Images= {'hist1.jpg','hist2.jpg','hist3.jpg'}
neg2cero=false;
for i=1:length(Images)
    I=imread(Images{i});
    [m,n,c]=size(I);
    [M,C] = MB_BKSVD4SD(double(I),D0,ns);
    C_Rmax= prctile(C',99);
    norm_fac=Cref_Rmax./C_Rmax;
    [Irec_norm,~]=Normaliza(Mref,C,M,neg2cero,norm_fac);
    Irec_norm=uint8(Irec_norm);
    subplot(2,3,i),imshow(I)
    subplot(2,3,i+3),imshow(Irec_norm)
    
    %Save results for posterior use
     patch_name = erase(Images{i},'.jpg');
     SaveResults(C,M,m,n,'./', patch_name)   

end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FOR PREVIOUSLY PATCHES IMAGES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Reference image
%Using patches
path_imgs='./hist1/*.jpg';
patches=dir(path_imgs);
folders={patches.folder};
patches={patches.name};
path_imgs=patches;
patches=cellfun(@(c)[folders{1} filesep c],patches,'uni',false);

% Estimates using all patches
[Mref,Cref_Rmax] = MI_MB_BKSVD4SD(patches,D0,ns);


%% image to normalize
%Using patches
path_imgs='./hist2/*.jpg';
patches=dir(path_imgs);
folders={patches.folder};
patches={patches.name};
path_imgs=patches;
patches=cellfun(@(c)[folders{1} filesep c],patches,'uni',false);

% Estimates using all patches
[M,C_Rmax] = MI_MB_BKSVD4SD(patches,D0,ns);
%% Deconvolución + normalizacion por parches

norm_fac=Cref_Rmax./C_Rmax;
for i=1:5
    loaded_images{i}=imread(patches{i});
end

figure(6)
neg2cero=false;
for i=1:length(loaded_images)
    
    CT=directDeconvolve(double(loaded_images{i}),M);
    [Irec_norm,~]=Normaliza(Mref,CT,M,neg2cero,norm_fac);
    Irec_norm=uint8(Irec_norm);

    subplot(2,length(loaded_images),i),imshow(loaded_images{i})
    subplot(2,length(loaded_images),i+length(loaded_images)),imshow(Irec_norm)
end

