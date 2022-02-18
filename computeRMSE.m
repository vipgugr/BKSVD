%
% Compute the RMSE (Elad's way)
% should split blocks of data to make computation less memory intensive
%

function err = computeRMSE(D,X,Y)


err2 = zeros(1,size(Y,2));
blocksize = 2000;
for i = 1:blocksize:size(Y,2)
  blockids = i : min(i+blocksize-1,size(Y,2));
  err2(blockids) = sum((Y(:,blockids) - D*X(:,blockids)).^2);
end

err = sqrt(sum(err2)/numel(Y));

end


