function [qtrain,qtest,softm_train,softm_test,alpha,D,U,R,C,affinity,V] = ICD(C,affinity,pivot,numclusters,Xtrain,Xtest,sim_type,sim_par)

%Clustering via the Incomplete Cholesky Decomposition with a stopping criterion based on 
%the stabilization of the cluster structure. The number of clusters is
%given

%Authors: Rocco Langone and Marc Van Barel

D = affinity;
D = D + eps;
D = 1./sqrt(D);
C2 = C;
C2 = bsxfun(@times,D,C2);
[Q,R] = qr(C2,0);
RR = R*R';
[V,E] = eig(RR);
diagE = diag(E);
[~,I] = sort(diagE);
V = V(:,I);
alpha = Q*V(:,end-numclusters+1:end);
U = V(:,end-numclusters+1:end);

%% Compute cluster memberships for training points
[qtrain,softm_train,QV] = compute_ci(alpha,D);

%% Compute memberships for test points, if any
if (exist('Xtest','var') && ~isempty(Xtest))
    
    Xpivots = Xtrain(pivot,:); 
    [qtest,softm_test] = compute_OoS(R,QV,U,Xpivots,Xtest,sim_type,sim_par); 
    
else
    
    qtest = [];
    softm_test = [];
end

