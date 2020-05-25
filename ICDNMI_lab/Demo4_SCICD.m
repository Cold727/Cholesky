%Demo studying the cluster performance of the fast spectral clustering 
%algorithm based on the Incomplete Cholesky Decomposition

%Author: Rocco Langone

%Citation: 
%Rocco Langone, Marc Van Barel and Johan A. K. Suykens, Entropy-Based Incomplete Cholesky Decomposition 
%for a Scalable Spectral Clustering Algorithm: Computational Studies and Sensitivity Analysis
%Entropy, Special Issue on Information Theoretic learning, June 2016.

clear
clc
close all

addpath(genpath('util'))

%% Define similarity type
datastruct.sim_type = 'rbf_sim';

%% Generate or load dataset
%R15
load R15.txt
X = R15(:,1:end-1);
X = preprocess_ICD(X);
labels = R15(:,end);

%% Settings
sigma = 0.25;
numclusters = 15;
datastruct.sim_par = sigma;
N = size(X,1);
datastruct.Xtrain = X;
THR_stop = 10^-6; %convergence threshold

%% Run algorithm
tStart = tic;
%Run algorithm with the given number of clusters
[qtrain,pivots,smothed_diff,h,softm_train,qtest,softm_test]=ICD2(datastruct,numclusters,THR_stop);
time = toc(tStart);
Xpivots = X(pivots,:); %coordinates of the pivot elements

%% Evaluate performance
ARI = adjrandindex(qtrain,labels);
num_pivots = length(pivots);

%% Plot results
figure
scatter(X(:,1),X(:,2),25,qtrain);
hold on
plot(Xpivots(:,1),Xpivots(:,2),'r*');
xlabel('x_1');
ylabel('x_2');
title('Clustered data + pivots');
box on
grid on
