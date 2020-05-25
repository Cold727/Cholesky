function [I,smtest] = compute_OoS(R,QV,U,Xtest,Xpivots,sim_type,sigma)

%Compute cluster indicators for test points

%Authors: Rocco Langone and Marc van Barel

K = sim_matrix(Xpivots,sim_type,sigma,Xtest);
Qext = K/R; %K*inv(R);
Vext = Qext*U;
Sext = Vext*QV;
Sext = Sext';
numclusters = size(Sext,1);
min_vett = min(Sext);
sm = Sext + repmat(2*abs(min_vett),numclusters,1); 
smtest = sm./repmat(sum(sm),numclusters,1);%soft memberships
[~,I] = max(abs(Sext));
        
