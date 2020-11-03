function  [CupDepth] = EvaluateCupDepth
   
   if ispc
      savepath = '../Save/Cup_2901/';
   else
      savepath = '/shared/OA/QuantifyCartilage/Save/Smoothness/';
   end
   
%    [Knees, KL1, KL2] = GetKL(CCBR);
   load('C:\Documents and Settings\stu.SUDHAKARTUMMALA\Desktop\Congruity\Matlab\ProgressionAnalysis\overviewFull.mat');
   
   [CupDepth] = LoadCurvatureFlowSmoothness(savepath, Knees);

   ev = Evaluate;
   
   Precision(CupDepth, 'Precision')
   separateOA(ev, CupDepth(:,1), KL(:,1), 'Cup Depth', 1, 2, 0, 1);

function [CupDepth] = LoadCurvatureFlowSmoothness(savepath, Knees)
  % Loads all pre-computed curvature values. For each knee there are three potential visits, so three values.
  % For non-found visits, the value is set to NaN.
  if ~exist(savepath,'dir')
     error(['Could not find ',savepath])
  end
  mats = dir([savepath,'/*.mat']);
  count = length(mats);
  if ~count, pathstr, error(['No smoothness in ',savepath]), end
  kneeCount = length(Knees);
  CupDepth = zeros(kneeCount,3);
  CupDepth(:) = NaN;
  for m = 1:count
     file = mats(m).name;
     knee = [file(1:10),' ',file(30)]; % Becomes CCBR number and then L or R
     visit = str2double(file(end-4)); % Visit is 1/3 for baseline
     if visit==6, visit=2; end        % Visit is 6 for follow-up, we put that in the second place
     load([savepath,'/',file]);
     % Check if the knees has already been found
     idx = strmatch(knee,Knees,'exact');
     if isempty(idx)
        error(['Did not find knee: ',knee])
     end
     % Record smoothness info
%      disp(sprintf('Flow Smoothness %03d/%03d: %s (%f,%f,%f,%f)', m, count, knee, smoothnessCurvature, smoothnessCurvatureROIS,smoothnessCurvatureROIM,smoothnessCurvatureROIL))
     CupDepth(idx,visit) = cupdepth;
  end

function Precision(Val,tag)
  V1 = Val(:,1);
  V2 = Val(:,3);
  ok = ~isnan(V1) & ~isnan(V2);
  V1 = V1(ok);
  V2 = V2(ok);
  stds = std([V1'; V2']);
  means = mean([V1'; V2']);
  cvs = stds ./ means;
  fprintf('%s (on %d scans): RMS CV %4.1f%% Mean CV %4.1f%% \n',tag,length(means),100*sqrt(nanmean(cvs.^2)),100*nanmean(cvs));
