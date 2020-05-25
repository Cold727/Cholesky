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

%% Load dataset
%Three  overlapping Gaussians
N = 10^3;
[X,labels] = threeclusters(N); 
X = preprocess_ICD(X); %normalize and remove outliers

%% Settings
THR_stop = 10^-6; %convergence threshold
sigma = mean(selectbandwidth(X,'method','SROTD'))*size(X,2); %Silverman's rule
maxk = 10;
datastruct.Xtrain = X;
datastruct.sim_par = sigma;
datastruct.Xtest= [-6 10;4 -10];

%% Run algorithm
tStart = tic;
%Select number of clusters based on the eigengap heuristics
[C,affinity,pivots,sim_type,sim_par,numclusters]=sel_clu_ICD(datastruct,maxk,THR_stop);
%Run algorithm with the selected number of clusters
[qtrain,qtest,softm_train,softm_test,alpha,D,U,R,C,affinity,V]=ICD(C,affinity,pivots,numclusters,X,[],sim_type,sim_par);
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

figure
numclu = length(unique(qtrain));
for i=1:numclu
    
    subplot(numclu,1,i);    
    scatter(X(:,1),X(:,2),30,softm_train(i,:),'filled')
    xlabel('x_1');
    ylabel('x_2');
    title(['Probability to belong to cluster ',num2str(i)]);
    box on
    grid on
end
