%Demo showing how to extract a subset representative of the cluster
%structure based on the Incomplete Cholesky Decomposition

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
sim_type = 'rbf_sim';

%% Generate or load dataset
%Two spirals
N = 10^3;
[X,labels] = twospirals(N); 
X = preprocess_ICD(X);
N = size(X,1);

%% Settings
sigma = 0.05;
m = 100;
idx_subset = subset_ICD(X,sim_type,sigma,m);
X_subset = X(idx_subset,:);

%% Plot results
figure
subplot(2,1,1)
scatter(X(:,1),X(:,2),25,labels,'fill');
xlabel('x_1');
ylabel('x_2');
title('Original data');
box on
grid on
subplot(2,1,2)
scatter(X_subset(:,1),X_subset(:,2),25,labels(idx_subset,:),'fill');
xlabel('x_1');
ylabel('x_2');
title('Extracted Subset');
box on
grid on
