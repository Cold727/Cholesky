function omega =sim_matrix(Xtrain,sim_type, sim_par,Xt)

% Construct the similarity matrix

% Author: Rocco Langone

switch lower(sim_type)
    
    case 'rbf_sim'
        
        if nargin<4
           omega = exp(-sqdist(Xtrain', Xtrain')/(2*sim_par^2));
        else
           omega = exp(-sqdist(Xtrain', Xt')/(2*sim_par^2));
        end
         
        
    case 'cosinerbf_sim'
             
        if nargin<4
          omega=exp(-squareform(0.5*pdist(Xtrain,'cosine').^2)/sim_par(1));

        else
          omega=exp(-(0.5*pdist2(Xt,Xtrain,'cosine').^2)/sim_par(1));

        end;
        
    case 'chisquared_sim'
        
        if nargin<4,       
            omega=exp(-squareform(pdist(Xtrain,@chi2distance))/sim_par(1));
        else

            omega=exp(-pdist2(Xt,Xtrain,@chi2distance)/sim_par(1));
        end;
    
    case 'cosine_sim'
        
        if nargin<4
            omega = 1- squareform(pdist(Xtrain,'cosine')); 
            omega(isnan(omega)) = eps;
        else
            omega = 1- pdist2(Xt,Xtrain,'cosine');
            omega(isnan(omega)) = eps;
        end;
            
    
    case 'corrrbf_sim'
        
        if nargin<4
          omega=exp(-squareform(0.5*pdist(Xtrain,'correlation'))/sim_par(1));

        else
          omega=exp(-(0.5*pdist2(Xt,Xtrain,'correlation'))/sim_par(1));   
     
        end;  
        
     case 'linear_sim'
        
        if nargin<4
            omega = Xtrain*Xtrain';
        else
            omega = Xtrain*Xt';
        end;      
end


