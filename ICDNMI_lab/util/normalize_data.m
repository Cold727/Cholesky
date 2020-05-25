function out_n = normalize_data(M)

N = size(M,1);
m = size(M,2);

for i =1:m
    
    media = mean(M(:,i));
    sig = std(M(:,i));
    
    M(:,i) = (M(:,i)-media)/sig;
    
end

out_n = M;