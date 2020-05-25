function [C,affinity,pivot,sim_type,sim_par,numclu_curr,smoothed_diff,cut_curr,h,h2]=sel_clu_ICD(datastruct,numclusters,threshold_stop)

%Clustering via the Incomplete Cholesky Decomposition with a stopping criterion based on 
%the stabilization of the cluster structure. The number of clusters is selected by means
%of the eigengap heuristics

%Authors: Rocco Langone and Marc van Barel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Xtrain = datastruct.Xtrain;
try
    sim_par = datastruct.sim_par;
catch
    sim_par = [];
end

try
sim_type = datastruct.sim_type;
catch
sim_type = [];    
end

N = size(Xtrain,1);
diagM = ones(N,1);
affinity = zeros(N,1);
pivot = zeros(1,N);
diff_vett = zeros(N-1,1);
GNcut = zeros(numclusters,N);
smoothed_diff = diff_vett;

maxGc = 200;
C = zeros(N,maxGc);
ok = false;
conv_w = 4; %convergence window
Iprev = [];

for l=2:numclusters
    numclu_curr = l;
    Iprev_curr = repmat((1:numclu_curr)',round(N/numclu_curr),1);
    if(length(Iprev_curr)< N)
        Iprev_curr = [Iprev_curr;ones(N-length(Iprev_curr),1)];
    else
        Iprev_curr = Iprev_curr(1:N);
    end
    Iprev = [Iprev;Iprev_curr'];
end

i = 1;
istart = 0;
numclu_prev = 2;
THR_aff = threshold_stop;

max_iter = 100;
max_iter2 = maxGc - max_iter;
break_flag = 0;

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
       
    %Only solve eigenvalue problem when the approximated affinity is good enough
    if((min(affinity)/max(affinity) > THR_aff) || (i>=max_iter))
    
        istart = istart + 1;
        
        %% Eigenvalue problem
        D = affinity;
        D = D + eps;
        D = 1./sqrt(D);
        C2 = C;
      
        C2 = bsxfun(@times,D,C2);
        [Q,R] = qr(C2,0);
        RR = R*R';
        
        try
            [V,E] = eig(RR);
            diagE = diag(E);
            [Y,I] = sort(diagE);
            GNcut(:,i) = Y(end-numclusters+1:end);
            V = V(:,I);
            alpha = Q*V(:,end-numclusters+1:end);
            [I,~] = compute_ci_all(alpha,D);
            cut_curr = GNcut(end-1:-1:1,i);
            eigengaps = abs(diff(cut_curr));
            [~,numclu_curr] = max(eigengaps); 
            numclu_curr = max(numclu_curr + 1,2);
            
        catch ME
            
            I = Iprev;
            ok = false;
        end
        
        %% Calculate stopping criterion
        
        if(i>1)
            
        try
        diff_curr = nmi(Iprev(numclu_prev-1,:),I(numclu_curr-1,:));
        diff_vett(i,:) = diff_curr;
        catch
        diff_vett(i,:) = diff_vett(i-1,:); 
        end
           
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
        numclu_prev = numclu_curr;
         
        else
            
        Iprev = I;
        ok = false;    
        numclu_prev = numclu_curr;
        
        end
        
    end
    
    i = i + 1;
    
 end
    
end
C = C(:,1:i-1);
pivot = pivot(1:i-1);
smoothed_diff = smoothed_diff(i-istart:i-1,:);

%% Plot
plot_flag=0;
%plot_flag=1;
if(~isempty(smoothed_diff) && (plot_flag==1))
%Convergence
h = figure;
plot(i-istart:i-1,smoothed_diff,'-bs','LineWidth',2);
max_diff = max(smoothed_diff);
min_diff = min(smoothed_diff);
axis([i-istart i-1 min_diff max_diff])
hold on
plot(i-istart:i-1,1-threshold_stop*ones(istart,1),'g--','LineWidth',2);
xlabel('x');
ylabel('y');
title('t')
box on

%Model selection plot
h2 = figure;
plot(1:numclusters,[1;cut_curr(end:-1:1,:)],'ms--','LineWidth',2,'MarkerFaceColor','m');
xlabel('x');
ylabel('y');
title('t')

else
h = [];
h2 = [];
end