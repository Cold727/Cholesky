function [X I]=twospirals(N)


N=N/2;
theta = linspace(pi/2,4*pi+pi/2,N);

a=0;
b=1;
factor=1;


noise=factor*(rand(size(theta))-0.5);
r=a+b*theta+noise;


[x,y]=pol2cart(theta,r);

X=[x' y'];
X=[X;-X];

I = [ones(N,1) ; 2*ones(N,1)];



