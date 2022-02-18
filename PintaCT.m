function [] = PintaCT(CT,M)
save=true;

%% EXPERIMENTAL RGB RECONSTRUCTION WITH FUNCTIONS
D=M;
X=CT;
[ns,p]=size(X);
c=3;
m=sqrt(p);
n=m;

Hhat=D(:,1)*X(1,:);
Ehat=D(:,2)*X(2,:);


%Return to rgb domain and image shape
Hhat = reshape(od2rgb(Hhat)',m,n,c);
Ehat = reshape(od2rgb(Ehat)',m,n,c);
% StainsRGB{1}=Hhat;
% StainsRGB{2}=Ehat;

figure()
subplot(1,ns,1);
%imshow(uint8(Hhat(1251:1500,1251:1500,:)));
imshow(uint8(Hhat));
title('Hematoxilina')
subplot(1,ns,2);
%imshow(uint8(Ehat(1251:1500,1251:1500,:)));
imshow(uint8(Ehat));
title('Eosina')

if ns>2
    Bhat=D(:,3)*X(3,:);
    Bhat = reshape(od2rgb(Bhat)',m,n,c);
%     StainsRGB{3}=Bhat;
    subplot(1,ns,3);
    imshow(uint8(Bhat));
    

    title('B')

end

if save==true
    imwrite(uint8(Hhat),'H.png','BitDepth',16);
    imwrite(uint8(Ehat),'E.png','BitDepth',16);
    if ns>2
        imwrite(uint8(Bhat),'B.png','BitDepth',16);
    end
end
end
