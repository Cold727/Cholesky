function [qtrain,pivot,smoothed_diff,h,softm_train,qtest,softm_test]=ICD2(datastruct,numclusters,threshold_stop)

%Fast Spectral Clustering using the Incomplete Cholesky Decomposition with
%a stopping criterion based on the stabilization of the cluster structure

%Authors: Rocco Langone and Marc van Barel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Xtrain = datastruct.Xtrain;
try
    sim_par = datastruct.sim_par;
catch
    sim_par = [];
end

try
    Xtest = datastruct.Xtest;
catch
    Xtest = [];
end

sim_type = datastruct.sim_type;

N = size(Xtrain,1);
diagM = ones(N,1);
affinity = zeros(N,1);
Iprev = repmat((1:numclusters)',round(N/numclusters),1);
if(length(Iprev)< N)
    Iprev = [Iprev;ones(N-length(Iprev),1)];
else
    Iprev = Iprev(1:N);
end
pivot = zeros(1,N);
diff_vett = zeros(N-1,1);
smoothed_diff = diff_vett;

maxGc = 200;
C = zeros(N,maxGc);
ok = false;
i = 1;
conv_w = 4; %convergence window
THR_aff = threshold_stop;
istart = 0;
break_flag = 0;

max_iter = 100;
max_iter2 = maxGc - max_iter;

while ( i <= N) & ~ok
 
 %% Select current pivot
    
 if(i==1) 
     
    [~,Ip] = max(diagM);
    pivot(i) = Ip(1);
    ctest = sim_matrix(Xtrain(:,:),sim_type,sim_par,Xtrain(Ip(1),:));
    break_flag = 0;
    
 else
    
    % Avoid the selection of the same pivot two times  
     
    [Y,Ip] = sort(diagM,'descend');
    out_flag = 1;
    idx_curr = 1;
    
    while(out_flag)
        
        ctest = sim_matrix(Xtrain(:,:),sim_type,sim_par,Xtrain(Ip(idx_curr),:));
        b=1; 
            
        if(sum(pivot(1:i)==Ip(idx_curr))==0) 
            out_flag = 0;
            Y = Y(idx_curr);
            Ip = Ip(idx_curr);
            pivot(i) = Ip;
        else
            affinity(Ip(idx_curr)) = max(affinity);
            idx_curr = idx_curr + 1;
        end
       
        if(idx_curr > length(ctest))
           out_flag = 0;
           break_flag = 1;
        end
       
    end

 end
 
    if(break_flag==1)
        
        break;
        
    else
        
    for ii = 1:i-1
        ctest = ctest - C(:,ii)*C(Ip(1),ii);
    end
    
    ctest=ctest + eps;
    c = ctest./sqrt(ctest(Ip(1)));
    affinity = affinity + c*(c'*ones(N,1)); %efficient version to compute the diagonal of C * C'
    C(:,i) = c;
    diagM = diagM - c.^2;
        
    if((min(affinity)/max(affinity) > THR_aff) || (i>=max_iter))
     
        istart = istart + 1;
        
        %% Eigenvalue problem
        D = affinity;
        D = D + eps;
        D = 1./sqrt(D);
        C2 = C;
              
        C2 = bsxfun(@times,D,C2);
        [Q,R] = qr(C2,0);
        clear C2
        RR = R*R';
        
        try
            [V,E] = eig(RR);
            diagE = diag(E);
            [~,I] = sort(diagE);
            V = V(:,I);
            alpha = Q*V(:,end-numclusters+1:end);
            I = compute_ci_new2(alpha,D);
            clear V
            clear alpha
            clear E
            clear diagE
            
        catch 
            I = Iprev;
            ok = false;
        end
          
        %% Compute stopping criterion 
        diff_curr = nmi(Iprev,I);
        diff_vett(i,:) = diff_curr;
                  
        if(i<=conv_w)
           idx_curr = 1:i;    
        else
           idx_curr = i-conv_w+1:i;
        end
                       
        smoothed_diff(1:i,:) = smooth(diff_vett(1:i,:),conv_w);
        avg = mean(smoothed_diff(idx_curr,:));
        ok = 1-avg < threshold_stop;
        ok = ok && (i>=numclusters);
        ok = ok || (istart>=max_iter2); 
        
        Iprev = I;
                
    end
    
    i = i+1;
    
    end
    
end
%% Compute cluster memberships for training points
smoothed_diff = smoothed_diff(1:i-1,:);
pivot = pivot(1:i-1);
C= C(:,1:i-1);
D = affinity + eps;
D = 1./sqrt(D);
C2 = C;
C2 = bsxfun(@times,D,C2);
[Q,R] = qr(C2,0);
RR = R*R';
[V,E] = eig(RR);
diagE = diag(E);
[~,I] = sort(diagE);
V = V(:,I);
U = V(:,end-numclusters+1:end);
alpha = Q*U;
[qtrain,softm_train,QV] = compute_ci(alpha,D);

%% Compute memberships for test points, if any
if (exist('Xtest','var') && ~isempty(Xtest))
    
    Xpivots = Xtrain(pivot,:); 
    [qtest,softm_test] = compute_OoS(R,QV,U,Xpivots,Xtest,sim_type,sim_par); 
    
else
    
    qtest = [];
    softm_test = [];
end

%% Plot convergence
%plot_flag = 1;
plot_flag = 0;
if(~isempty(smoothed_diff(1:i-1,:)) && (plot_flag==1))
h = figure;
plot(1:i-1,smoothed_diff,'-bs','LineWidth',2);
max_diff = max(smoothed_diff);
min_diff = min(smoothed_diff);
axis([1 i-1 min_diff max_diff])
hold on
plot(1:i-1,1-threshold_stop*ones(length(smoothed_diff),1),'g--','LineWidth',2);
xlabel('x');
ylabel('y');
title('t')
box on
grid on
else
h = [];
end    
