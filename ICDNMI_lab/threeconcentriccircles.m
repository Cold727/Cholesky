function [X Ioriginal] = threeconcentriccircles(N)

N1 = floor(N/2);
N2 = floor(N/4);
N3 = N-N1-N2;
sumN = N1 + N2 + N3;
X = zeros(sumN,2);
factor = 0.2;
%circle 1
r1 = 4;%radius 
theta1 =linspace(0,2*pi,N1);
noise1 = factor*(rand(size(theta1)));
X(1:N1,1) = r1*cos(theta1)+noise1;
X(1:N1,2) = r1*sin(theta1)+noise1;

%circle 2
r2 = 1.8;
theta2 =linspace(0,2*pi,N2); 
noise2 = factor*(rand(size(theta2)));
X(N1+1:N1+N2,1) = r2*cos(theta2)+noise2;
X(N1+1:N1+N2,2) = r2*sin(theta2)+noise2;

%circle 3
r3 = 1.0;
theta3 = linspace(0,2*pi,N3);
noise3 = factor*(rand(size(theta3)));
X(N1+N2+1:sumN,1) = r3*cos(theta3)+noise3;
X(N1+N2+1:sumN,2) = r3*sin(theta3)+noise3;

Ioriginal = [ones(N1,1); 2*ones(N2,1); 3*ones(N3,1)];




