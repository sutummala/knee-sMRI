%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function N = SampleSized(Eval, H, D, alpha, varequal, verbose, desPWR)
% Uses Satterthwaite's approximation for the effective degrees of freedom
% for unequal variance
% When vartype is 2 (unequal), the test assumes that the two samples come from normal distributions
% with unknown and unequal variances. This is known as the Behrens-Fisher-Welch problem

if ~exist('alpha','var')
   alpha = 0.05;
end
if ~exist('varequal','var')
   varequal = 'equal';
end
if ~exist('verbose','var')
   verbose = 0;
end
if ~exist('desPWR','var')
   desPWR = 0.8;
end

N = SampleSize(H, D, alpha, varequal, verbose, desPWR);
