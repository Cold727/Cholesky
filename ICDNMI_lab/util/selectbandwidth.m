function h = selectbandwidth(varargin)

% Description
% This function selects the bandwidth of the RBF kernel for density
% estimation. The bandwidth can bve selected according to reference rules
% or a plug-in estimate.  
%
% The refence rules that can be chosen are:
% - SNR or snr     : Simple Normal Reference
% - SROT or srot   : Silverman's Rule of Thumb for one dimensional data
% - OS or os       : OverSmoothed bandwidth
% - SROTD or srotd : Silverman's Rule of Thumb for D dimensional data
% - ISNR or isnr   : Improved Simple Normal Reference based on Edgeworth
%                    expansion for f around the Gaussion density
%
% 'SNR' is the default for 1 dimensional data and 'SROTD' is the default
% for multi dimensional data.
%
% The plug-in estimator
% - SJPI or sjpi   : Sheather-Jones plug-in calculated via density
%                    functionals
%
% Full syntax
%
% >> h = selectbandwidth(X)
% >> h = selectbandwidth(X,'method',type)
%
% only when using 'sjpi'
% >> h = selectbandwidth(X,'method',type,'tol',value) 
%
% Outputs
%   h       : 1 x D bandwidth vector 
%
% Inputs
%   X       : N x D matrix with input data
%   type(*) : 'SNR(*)','SROT','OS','SROTD(*)' or 'SJPI'

nvar = numel(varargin);
X = varargin{1,1};
[n,d] = size(X);

if nvar <= 2 && d == 1 % default for 1D data
    opt = 'SNR';
elseif nvar <= 2 && d > 1 % default for multi D data
    opt = 'SROTD';
else
    opt = varargin{1,3};
end

switch opt
    
    case {'SNR','snr'}
        s = min(std(X),iqr(X)/1.349);
        h = 1.06*s*n^(-1/5);
    
    case {'SROT','srot'}
        A = min(std(X),iqr(X)/1.349);
        h = 0.9*A*n^(-1/5);
        
    case {'OS','os'}
        s = std(X);
        h = 1.144*s*n^(-1/5);
        
    case {'SROTD','srotd'}
        s = std(X);
        h =  s*(4/(n*d+2*n))^(1/(d+4));
        
    case {'SJPI','sjpi'}
        if nvar < 4
            tol = 1e-3;
        else
            tol = varargin{1,5};
        end
        h = zeros(1,d);
        for i=1:d
            h(i) = selectbandwidth(X(:,i),tol);
        end
        
    case{'ISNR','isnr'}
        c1=35/48; c2=35/32; c3=385/1024;
        g3=skewness(X,0); g4=kurtosis(X,0);
        s = min(std(X),iqr(X)/1.349);
        h_snr = 1.06*s*n^(-1/5);
        h = h_snr.*(1+c1*g4+c2*g3.^2+c3*g4.^2).^(-1/5);

    otherwise
        error('Unknown method. Please choose ''SNR'',''SROT'',''SROTD'',''OS'', ''SJPI'' or ''ISNR''.')
end


 