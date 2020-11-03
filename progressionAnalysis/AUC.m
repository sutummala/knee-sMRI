function a = AUC(scr,lab,visualize)

if ~exist('visualize','var'), visualize=0; end

[uscr,a,map] = unique(-scr); % Unique sorts in the wrong order, so "-" handles that

% [scr,idx] = sort(scr);
% lab = lab(idx);
len = numel(uscr);
ROC = zeros(len+2,2);
pos = sum(lab==1);
fal = sum(lab==0);
tp = 0;
fp = 0;
for s = 1:len
   labidx = find(map==s);
   labs = lab(labidx);
   tp = tp+sum(labs==1);
   fp = fp+sum(labs==0);
   ROC(s+1,:) = [fp/fal,tp/pos];
end

% Add points for thresholds -inf and +inf
if ROC(2,2)<ROC(end-1,2)
   ROC( 1 ,:) = [0 0];
   ROC(end,:) = [1 1];
else
   ROC(end,:) = [0 0];
   ROC( 1 ,:) = [1 1];
end

if visualize
   % disp(ROC')
   figure, plot(ROC(:,1),ROC(:,2),'r-*')
   ylabel('True Positive')
   xlabel('False Positive')
end

% Compute AUC
ROC(end+1,:) = [1,0];
a = polyarea(ROC(:,1),ROC(:,2));

if a < 0.5
   a = 1-a;
end
