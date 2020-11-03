function [A, Cov] = DeLongCovarians(X,Y)

% X include biomarkervalues for normals 
% Y include biomarkervalues for abnormals

[no_markers no_normals] = size(X); 
[no_markers no_abnormals] = size(Y); 

% Compute structural components

vab = zeros(no_markers,no_abnormals);

for m = 1:no_markers
    for i = 1:no_abnormals
        vab(m,i) = (sum(X(m,:) < Y(m,i)) + 0.5 * sum(X(m,:) == Y(m,i)))/no_normals;
    end;
end;


vno = zeros(no_markers,no_normals);

for m = 1:no_markers
    for i = 1:no_normals
        vno(m,i) = (sum(X(m,i) < Y(m,:)) + 0.5 * sum(X(m,i) == Y(m,:)))/no_abnormals;
    end;
end;

% Compute estimate of areas

A = zeros(no_markers,1);
for m = 1:no_markers
    A(m) = sum(vab(m,:))/no_abnormals;
end;

% compute the estimated covariance matrix of areas

for m = 1:no_markers
    for n = 1:no_markers
        s10m(m,n) = (vab(m,:)-A(m)) * (vab(n,:)-A(n))' / (no_abnormals-1);
      
    end;
end;

for m = 1:no_markers
    for n = 1:no_markers
        s01m(m,n) = (vno(m,:)-A(m)) * (vno(n,:)-A(n))' / (no_normals-1);
      
    end;
end;

Cov = s10m/no_abnormals + s01m/no_normals;

