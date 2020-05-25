function [X I]=threeclusters(n)

nsmall=floor(n/3);


c1=[-8,2.5];

rscaling = 0.5;

X=zeros(n,2);
for i=1:nsmall
    r=randn(1)/rscaling;
    theta = rand(1)*2*pi;   
    X(i,1) = c1(1)+r*sin(theta);
    X(i,2) = c1(2)+r*cos(theta);    
end

%c2 = [3,6.5];
c2 = [4,10];

for i=nsmall+1:2*nsmall
    r=randn(1)/rscaling;
    theta = rand(1)*2*pi;   
    X(i,1) = c2(1)+r*cos(theta);
    X(i,2) = c2(2)+r*sin(theta);
    
end   


c2 = [0,-10];
for i=2*nsmall+1:n
    r=randn(1)/rscaling;
    theta = rand(1)*2*pi;    
    X(i,1) = c2(1)+r*sin(theta);
    X(i,2) = c2(2)+r*cos(theta);    
end  
 I = [ones(1,nsmall) 2*ones(1,nsmall) 3*ones(1,n-2*nsmall)];

 %plot(X(:,1),X(:,2),'o')