
function [CT] = directDeconvolve(I,RM)

    Y=rgb2od(I);
    [m,n,c]=size(Y);
    YT=reshape(Y,(m)*(n),c)';
    CT = RM \ YT;
    
end