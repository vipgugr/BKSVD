function [] = PintaMatriz(M)
figure()
ns=size(M,2)
Mrgb=od2rgb(M);

for i=1:ns
    H=ones(3,10000);
    H=H.*Mrgb(:,i);
    H=reshape(H',[100,100,3]);
    subplot(1,ns,i)
    imshow(uint8(H))
end