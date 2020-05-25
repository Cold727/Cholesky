function [I, Y, QV, repres] = compute_ci_all(alpha,D)

%Computation of cluster indicators for a range of cluster numbers [2:k]

%Authors: Rocco Langone and Marc van Barel

numclusters = size(alpha,2);
DD = D*ones(1,numclusters);
alpha = alpha./DD;
alpha = alpha(:,end:-1:1);
alpha = alpha';
I = zeros(numclusters-1,size(alpha,2));
Y = zeros(numclusters-1,size(alpha,2));

for k = 2:numclusters
    [QV,RV,EV]=qr(alpha(1:k,:),0);
    RR = RV(1:k,1:k)\RV;
    RR(:,EV) = RR;
    [Y(k-1,:),I(k-1,:)] = max(abs(RR));
    repres = EV(1:numclusters);
end

