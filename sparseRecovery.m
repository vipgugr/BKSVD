function [X, Y_rec] = sparseRecovery(D,Y,nmi)

Q = size(Y,2);

disp('Computing sparse representation:')

%pbar(0,'initialize')

X=zeros(size(D,2),Q);

if nargin == 2
    s = 0;
    for q=1:Q
        [weights,~,used,~,~,~,~,~,~] = FastLaplace_add(D,Y(:,q));
        X(used,q)=weights;
        s = pbar2(q/Q,s);
    end
    Y_rec = D*X;
else
    s = 0;
    for q=1:Q
        [weights,~,used,~,~,~,~,~,~] = FastLaplace_add(D(nmi(:,q),:),Y(nmi(:,q),q));
        X(used,q)=weights;
        s = pbar2(q/Q,s);
    end
    Y_rec = D*X;
end