%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function [OR, threshold] = OddsRatio(H, D, verbose)

  thresholds = 100;

  if ~exist('verbose','var')
     verbose = 1;
  end

  % To make it simple, we swap to Diseased are always elevated values
  if mean(H)>mean(D)
     temp = D;
     D = H;
     H = temp;
  end

  if 1 % Pick median value and compute OR
      threshold = (median([H,D]));
      OR = Ratio(H,D,threshold);
  else % Find optimal threshold (between means) and OR
      % We assume that the distributions are Gaussian-like. So the
      ts = linspace(mean(H),mean(D),thresholds);
      OR = -1;
      threshold = NaN;
      for t = 1:thresholds
          ort = Ratio(H,D,ts(t));
          if ort>OR
              OR = ort;
              threshold = ts(t);
          end
          orts(t) = ort;
      end
      if verbose
          figure, plot(ts,orts)
      end
  end

function OR = Ratio(H, D, threshold)
  % Count of values Higher than threshold
  HF = length(find(H>=threshold)); % false healthy
  DT = length(find(D>=threshold)); % true diseased
  % Count of values Lower than threshold
  HT = length(find(H<threshold));
  DF = length(find(D<threshold));
  % Odds Ratio
  if HF==0 && DF==0
     OR = Inf; % Perfect!
  elseif HT==0 || DT==0
     OR = 0; % Horrible!
  else
     OH = HT / HF;
     OD = DT / DF;
     OR = OH * OD;
     if 0 & OR < 1
        disp('OddsRatio: something weird here?!')
        OR = 1 / OR;
        % Because we dont know if the threshold is upper or lower.
        % And since we only use OR for groups with significant, Gaussian
        % non-overlap, then there will be a correct OR over 1.
     end
  end
