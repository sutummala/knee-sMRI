function [orts, stars] = RatioVsLow(H, D, quantWidth, makeFig, offset)
  global useColor
  if ~exist('makeFig','var'), makeFig=0; end
  if ~exist('offset' ,'var'), offset=[]; end

  quantLows = 0.0:quantWidth:0.98;
  thresholds = numel(quantLows);
  ts = 0*quantLows;
  for q = 1:thresholds
     tsLo(q) = quantile([H D],quantLows(q));
  end
  tsLo(1)=-Inf; tsLo(end+1)=Inf;
  if quantWidth==0.5 && numel(unique([H D]))==2
     % Hack to ensure binary data does not confuse the equalities below
     tsLo(2) = mean([H D]);
  end

  % To make it simple, we swap to Diseased are always elevated values
  if mean(H)>mean(D)
     temp = D;
     D = H;
     H = temp;
  end
  
  OR = -1;
  orts = 0*thresholds;
  stars = orts; CIlo = orts; CIhi = orts;
  for t = 1:thresholds
     % Standard Odds Ratio
     [ort,MHOR,MHCI,star] = Ratio(H,D,tsLo(2),tsLo(t),tsLo(t+1));
     orts(t) = ort;
     stars(t) = star;
     CIlo(t) = MHCI(1);
     CIhi(t) = MHCI(2);
  end

  if ~nargout
     fprintf('Ratios vs low: ');
     for t = 2:thresholds
        fprintf('%f %.1f ',tsLo(t),orts(t));
     end
     fprintf('\n');
  end
  
  if makeFig
     if isempty(offset) || offset==0;
        figure
        set(gca,'fontsize',20)
        hold on
     end
     xs = 100*(quantLows+0.5*quantWidth);
     if isempty(offset)
        h = bar(xs, orts);
        if useColor, col = [102 149 160] / 256;
        else         col = [0.6 0.6 0.6]; end
        offset = 0;
     else
        h = bar(xs+offset*100, orts, 0.3);
        if useColor, col = [80 130 140] / 256 +2*offset;
        else         col = [0.4 0.4 0.4]+2*offset; end
        end
     set(h,'facecolor',col)
     % errorbar(xs, orts, CIlo, CIhi, 'linestyle','none','linewidth',10,'color','k')
     ylabel('Odds Ratio')
     yl = ylim;
     ylim([0 yl(2)]);
     yl = 0.05 * yl(2);
     for t = 1:thresholds
        if stars(t)
           text(xs(t)+offset*100,orts(t)+yl,repmat('*',1,round(stars(t))),'horizontalalignment','center','fontsize',28)
        end
     end
     switch quantWidth
        case 0.25, set(gca,'xticklabel',{'Q1','Q2','Q3','Q4'}), xlabel('Quartiles')
        case 0.33, set(gca,'xtick',[16.6 50 83.3], 'xticklabel',{'T1','T2','T3'}), xlabel('Tertiles')
        case 0.5,  set(gca,'xtick',[25 75],'xticklabel',{'Low','High'})
        otherwise, xlabel('Fractile'),
     end
  end
  

function [OR,MHOR,MHCI,stars] = Ratio(H, D, bottom, threshLo, threshHi)
  % Count of values inside thresholds
  HF = length(find(H>=threshLo & H<threshHi)); % false healthy
  DT = length(find(D>=threshLo & D<threshHi)); % true diseased
  % Count of values outside thresholds
  HT = length(find(H<bottom));
  DF = length(find(D<bottom));
  table = [HT HF; DF DT];
  OR = OddsRatio(table);
  [MHOR,MHCI,stars] = MHratio(table);
  
function OR = OddsRatio(table)
  HF = table(2);
  DT = table(4);
  HT = table(1);
  DF = table(3);
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

function [MHOR,MHCI,stars] = MHratio(table)
  alpha = 0.05;
  [MHOR varlogtheta MHCI] = MantelHaenszelOddsRatio(table,alpha);
  stars = 0;
  if MHCI(1)>1
     stars = 1; alpha = 0.01;
     [OR varlog CI] = MantelHaenszelOddsRatio(table,alpha);
     if CI(1)>1
        stars = 2; alpha = 0.001;
        [OR varlog CI] = MantelHaenszelOddsRatio(table,alpha);
        if CI(1)>1
           stars = 3; alpha = 0.0001;
           [OR varlog CI] = MantelHaenszelOddsRatio(table,alpha);
           if CI(1)>1
              stars = 4;
           end
        end
     end
  end
  
  
function [theta varlogtheta ci]  = MantelHaenszelOddsRatio(x,alpha)
% # Function to produce Mantel-Haenszel common odds ratio estimate and asymptotic CI
% # Function based on Thompson (1999)
% # Input : x - 2 x 2 x K array of counts
% # Output: theta - MH estimate
% #         ci - 100 x (1-alpha)% confidence interval for theta

n11k = x(1,1,:);
n22k = x(2,2,:);
n12k = x(1,2,:);
n21k = x(2,1,:);
nk = n11k+n22k+n12k+n21k;

theta = sum(n11k.*n22k./nk)/sum(n12k.*n21k./nk);
varlogtheta = sum( (n11k+n22k).*n11k.*n22k./(nk.*nk)) /(2*sum(n11k.*n22k./nk)*sum(n11k.*n22k./nk)) + ...
    sum(((n11k + n22k).*n12k.*n21k + (n12k + n21k).*n11k.*n22k)./(nk.*nk))/ (2*sum(n11k.*n22k./nk)*sum(n12k.*n21k./nk)) + ...
    sum((n21k+n21k).*n12k.*n21k./(nk.*nk)) / (2*sum(n12k.*n21k./nk)*sum(n12k.*n21k./nk));
zalpha = norminv(1-alpha/2,0,1);
ci = [exp(log(theta)-zalpha*sqrt(varlogtheta))   exp(log(theta) + zalpha*sqrt(varlogtheta))];
