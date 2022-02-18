function [Inorm, StainNorm] = Normaliza(MT,C,M,neg2cero,norm_fac)

if neg2cero
    C(C < eps) = eps;
end


% D=M;
% X=C;
[ns,p]=size(C);
c=3;
m=sqrt(p);
n=m;

CT_norm=C.*norm_fac';
Yrec=MT(:,1:ns)*CT_norm;
Y2d=reshape(Yrec',m,n,c);
%Irec=OD2intensities(Y2d);
Inorm=od2rgb(Y2d);

for i=1:ns
    stain=MT(:,i)*CT_norm(i,:);
    stain=reshape(od2rgb(stain)',m,n,c);
    StainNorm{i}=stain;
end



