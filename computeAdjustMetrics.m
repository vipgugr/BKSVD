%
% Compute the RMSE (Elad's way)
% should split blocks of data to make computation less memory intensive
%

function [BIC,AIC] = computeAdjustMetrics(D,X,Y)


err2 = zeros(1,size(Y,2));
blocksize = 2000;
for i = 1:blocksize:size(Y,2)
  blockids = i : min(i+blocksize-1,size(Y,2));
  err2(blockids) = sum((Y(:,blockids) - D*X(:,blockids)).^2);
end

%err = sqrt(sum(err2)/numel(Y));
serr2=sum(err2);
n=size(Y,2);
ns=size(X,1);
L=serr2/n;

%BIC=n*log(serr2/n)+ns*log(n);
%AIC=n*log(serr2/n)+2*ns;

% Cp=(1/n)*(serr2+2*ns);
% BIC=(1/n)*(serr2+log(n)*ns);
% AIC= (1/n)*(serr2+2*ns)

BIC=log(n)*ns+2*log(L);
AIC=2*ns-2*log(L);

end

%%
% Y=Y5;
% D=D2{5};X=X2{5};
% prueba=[69.88 1.42e-13; 40.27 9.71e-14; 27.85 3.5e-14]/3*n
% ML=[0.0049 3.87e-32; 0.0016 1.83e-32; 0.003 9.51e-32]/3
% BIC=zeros(3,2);
% BIC(:,1)=log(n)*2;
% BIC(:,2)=log(n)*3;
% BIC+2*log(ML)
% BIC+2*log(prueba)