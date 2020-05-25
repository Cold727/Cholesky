function Xnew = preprocess_ICD(X)

%Remove outliers and scale data 

%Author: Rocco Langone

Xtemp = normalize_data(X);
med_temp = median(Xtemp);
std_temp = std(Xtemp);
I = bsxfun(@gt, abs(bsxfun(@minus, Xtemp, med_temp)), 3*std_temp);
I = logical(sum(I,2));
Xtemp(I,:) = [];
Xnew = Xtemp;

