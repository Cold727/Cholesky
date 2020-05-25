function [I,sm,QV,repres] = compute_ci(alpha,D)

%Computation of cluster indicators for a given number of cluster k

%Authors: Rocco Langone and Marc van Barel

numclusters = size(alpha,2);

if(numclusters~=1)

    DD = D*ones(1,numclusters);
    alpha = alpha./DD;
    alpha = alpha(:,end:-1:1); 
    alpha = alpha';
    [QV,RV,EV]=qr(alpha,0);
    RR = RV(1:numclusters,1:numclusters)\RV;
    RR(:,EV) = RR; 
    min_vett = min(RR);
    sm = RR + repmat(2*abs(min_vett),numclusters,1);
    sm = sm./repmat(sum(sm),numclusters,1);
    [~,I] = max(abs(RR));
    repres = EV(1:numclusters);
    
else
    N = size(alpha,1);
    I = ones(1,N);
    sm = I;
    DD = D*ones(1,numclusters);
    alpha = alpha./DD;
    alpha = alpha(:,end:-1:1); 
    alpha = alpha';
    QV=qr(alpha,0);
    repres = [];
end