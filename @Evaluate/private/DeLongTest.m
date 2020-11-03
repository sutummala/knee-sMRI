function [dA, sddA, p] = DeLongTest(X,Y,L)
 
% X are the normals
% Y are the abnormals
% L is the contrast function 

% In case of two biomarkers
%    X(1,:) are marker 1 values
%    X(2,:) are marker 2 values
%    Y(1,:) are marker 1 values
%    Y(2,:) are marker 2 values

% testing Aroc(M1) > Aroc(M2) -> L = [1 -1];

% For 3 markers testing Aroc(M1) larger than any other
%      ->   L = [[1 -1 0]' [1 0 -1]']';

[A Cov] = DeLongCovarians(X,Y);

%A
dA = L'*A;

VardA = L' * Cov * L;
sddA=sqrt(VardA);

if (size(X,1) > 1)
  p = 1- chi2cdf(dA' * inv(VardA) * dA,rank(L));
else
  p = 2*normcdf(0.5, A, sddA);
  if (p>1)
    p = 2-p;
  end
end

