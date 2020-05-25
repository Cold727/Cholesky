function pivot=subset_ICD(Xtrain,kernel_type,sigmas,m)

%Extract a subset consisting of the Incomplete Cholesky Decomposition pivots, where the ICD 
%approximates the kernel matrix

%Author: Rocco Langone 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = size(Xtrain,1);
num_sigmas = length(sigmas);
pivot = zeros(m,num_sigmas);

for j=1:num_sigmas

    diagM = ones(N,1);
    C = zeros(N,m);
    
for i = 1:m
 
 %Select current pivot
 [~,Ip] = max(diagM);
 pivot(i,j) = Ip(1);
 ctest = sim_matrix(Xtrain,kernel_type,sigmas(j),Xtrain(Ip(1),:));
        
 for ii = 1:i-1
        ctest = ctest - C(:,ii)*C(Ip(1),ii);
 end
    
 ctest=ctest + eps;
 c = ctest./sqrt(ctest(Ip(1)));
 C(:,i) = c;
 diagM = diagM - c.^2;
  
end

end

[d1,d2] = size(pivot);
pivot = unique(reshape(pivot,d1*d2,1));
