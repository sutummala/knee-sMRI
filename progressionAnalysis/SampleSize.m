%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function N = SampleSize(H, D, alpha, varequal, verbose, desPWR)
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
if strcmp(varequal,'unequal')
   vartype = 2;
else
   vartype = 1;
end
dim = 1;
a = alpha;

x = H(:);
y = D(:);

nx = numel(x);
ny = numel(y);
% Variances
s2x = nanvar(x, [], dim);
s2y = nanvar(y, [], dim);
% Difference in the means
difference = nanmean(x, dim) - nanmean(y, dim);

if vartype == 1 % equal variances
    % Degrees of freedem for equal variance
    dfe = nx + ny - 2;
    % Pooled standard error
    sPooled = sqrt(((nx-1) .* s2x + (ny-1) .* s2y) ./ dfe);
    se = sPooled .* sqrt(1./nx + 1./ny);
    ratio = difference ./ se;

    stats = struct('tstat', ratio, 'df', cast(dfe,class(ratio)), 'sd', sPooled);
    if isscalar(dfe) && ~isscalar(ratio)
        stats.df = repmat(stats.df,size(ratio));
    end
elseif vartype == 2 % unequal variances
    s2xbar = s2x ./ nx;
    s2ybar = s2y ./ ny;
    % Degrees of freedom using Satterthwaite's approximation
    dfe = (s2xbar + s2ybar) .^2 ./ (s2xbar.^2 ./ (nx-1) + s2ybar.^2 ./ (ny-1));
    se = sqrt(s2xbar + s2ybar);
    % T-statistic or effect size, Behrens-Welch 
    ratio = difference ./ se;

    stats = struct('tstat', ratio, 'df', cast(dfe,class(ratio)), 'sd', sqrt(cat(dim, s2x, s2y)));
    if isscalar(dfe) && ~isscalar(ratio)
        stats.df = repmat(stats.df,size(ratio));
    end
end

% Comput a two tailed p-value, thats why multiply by 2
% tcdf Student's T cumulative distribution function, ratio is effect size
% and dfe is degrees of freedom
p = 2 * tcdf(-abs(ratio), dfe);

if verbose
   if vartype == 1
      disp(sprintf('p (equal variance) = %e', p));
   else
      disp(sprintf('p (unequal variance) = %e', p));
   end
end

a = a/2;
tb1 = ratio - tinv(1 - a, dfe);  %Power estimation.
tb2 = ratio + tinv(1 - a, dfe);
 b1 = 1 - tcdf(tb1, dfe);
 b2 = 1 - tcdf(tb2, dfe);
Power = 1 - (b1 - b2);
if verbose
   disp(sprintf('Power is: %2.4f', Power));
end

% va is the average variance of the two groups
va = (s2x + s2y)/2;
% This equation gives the number of N for alpha level and power, multiply
% by 4 to get for both sets else by 2 for one set, norminv is Inverse of the
% normal cumulative distribution function
N = ceil( 4*va*(norminv(1-a) + norminv(desPWR))^2 / difference^2 );
if verbose
   disp(sprintf('Total N for 0.8 power: %d', N));
end
