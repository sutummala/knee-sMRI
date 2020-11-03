function  [Congruity] = EvaluateCongruity
   
   if ispc
      savepath = '../Save/Measure_all_0512_1_manual/';
      %savepath = '../Save/Measure_all_1801_1/';
      % savepath = '../Save/training25/';
   else
      savepath = '/shared/OA/QuantifyCartilage/Save/congruity/';
   end
   
%    [Knees, KL1, KL2] = GetKL(CCBR);
   load('C:\Documents and Settings\stu.SUDHAKARTUMMALA\Desktop\Congruity\Matlab\ProgressionAnalysis\overviewFull.mat');
   
    removeKnees = []; 
    for i = 1:length(Knees)
        x = strcmp(Knees(i), sort(training25));
        if sum(x), removeKnees(end+1) = i; end;
    end
    
   ExcludedKL = find(isnan(SmoothT(:,1))); removeKnees = union(ExcludedKL, removeKnees); removeKnees = union(removeKnees,63);
   ok = setdiff(1:318, removeKnees);
   [Congruity] = LoadCongruity(savepath, Knees,2,2,'Manual');
   
   %% Remove the values for the knees used for training of the clsssifier
   Congruity = Congruity(ok,:); nVolTF = nVolTF(ok,:); nThickT = nThickT(ok,:); SmoothT = SmoothT(ok,:); % MRI markers
   Interval = Interval(ok); % Baseline - Followup interval
   KL = KL(ok,:); nGap = nGap(ok,:); %Radiograph markers
   Age = Age(ok,:); BMI = BMI(ok,:); %Confounding factors
   CTX2cr = CTX2cr(ok,:); % Cartilage degradation biochemical marker
   Sex = Sex(ok,:); men = Sex(:,1) == 0; women = Sex(:,1) == 1; % Gender/Sex
   Pain = Pain(ok,:); %Clinical
   %-----------------------------------------------------------
   
   %% Correct for Age and BMI
   
   correct = 0;
   if correct
      [Cong] = correctSomething(Congruity, KL, Age, BMI, Sex);
      separateOA(ev, Cong(:,1), KL(:,1), 'Congruity Index [mm]', 1, 1, 1, 1);
   end
   
   %% After Menopause, particularly interesting while looking at Gender differences
   mPause = 0;
   if mPause
        menopause = Age(:,1) >= 50;
        Congruity = Congruity(menopause,:); KL = KL(menopause,:); Sex = Sex(menopause,:);
        men = Sex(:,1) == 0; women = Sex(:,1) == 1;
   end
   ev = Evaluate;
   
   %% Age, BMI predicts/correlated with Congruity---------
   
   agebmi = 0;   
   if agebmi
   healthy = KL(women,1) >= 0;
   age = BMI(women,1); cong = Congruity(women,1);
   h = chi2gof(cong);
    if ~h
            [c1, p1] = corr(age(healthy), cong(healthy), 'type', 'Pearson');
    else
            [c2, p2] = corr(age(healthy), cong(healthy), 'type', 'Spearman');
    end
   end
   
   
   %% both baseline and follow-up for SRM---------------------------
   
   srm = 0;
   if srm
        basfol = ~isnan(Congruity(:,1))& ~isnan(Congruity(:,2));
        srm = (Congruity(basfol,2) - Congruity(basfol,1))./Interval(basfol); % This is not SRM. don't be panic!
        lChange = 100 * srm./(Congruity(basfol,1)./Interval(basfol));
        KL = KL(basfol,:);
        men = Sex(basfol,1) == 0; women = Sex(basfol,1) == 1;
   end
   %--------------------------------------------------------------
   
   %% With pain----------------------------------------
   pai = 0;
   if pai
    ok1 = KL(men,2) > 1; % For '0' healthy, '1' early OA and for '>1' OA
    pain = Pain(men,2); con = Congruity(men,2);
    subpain = pain(ok1); subcon = con(ok1);
    con1 = subpain(subpain ~= 0); con2 = subcon(subpain ~= 0);
    notnan = ~isnan(con1) & ~isnan(con2);
    [c,p] = corrcoef(con1(notnan), con2(notnan)); % Use for correlations
%     [OR,stars] = RatioVsLow(con1',con2',0.50);
%     [h,p1] = ttest2(con1, con2);
   end
   %---------------------------------------------------
   
   both = Precision(Congruity(:,:), 'Precision'); % precision
   
   %% Correaltions between volume loss/JSN and congruity change
   efficacy = 0;
   if efficacy
    okey = KL(men,1) == 3;
    edi = 2; cong = Congruity(okey,1);
    dif = longchange(nGap(men,:), edi, okey); % Volume Change
    dif1 = longchange(Congruity(men,:), edi, okey);% Congruity Change
   
    not = ~isnan(dif) & ~isnan(dif1) & ~isinf(dif);
    h1 = chi2gof(dif); h2 = chi2gof(dif1); % Checking whether the data is normally distributed
   
    if ~h1 && ~h2
            [cor, p] = corr(dif(not), dif1(not), 'type', 'Pearson');
    else
            [cor1, p1] = corr(dif(not), dif1(not), 'type', 'Spearman');
    end
   end
   
   %% Whether Congruity predicts cartilage loss
   predictOA = 0;
   if predictOA
    if 1
            prog = Congruity((dif < nanmedian(dif)),1);
            nonprog = Congruity((dif >= nanmedian(dif)),1);
    else
            prog = dif(cong < nanmedian(cong));
            nonprog = dif(cong >= nanmedian(cong));
    end
   
   perpoints = (nanmean(prog) - nanmean(nonprog)) * 100;
   [h,p] = ttest2(prog, nonprog);
   end
   
   %Graph showing the Cross-sectional separation based on the KL index
   separateOA(ev, Congruity(:,1), KL(:,1), 'Congruity Index', 1, 1, 1, 1); % Graph showing seperation of KLs
   
   
function [Congruity] = LoadCongruity(savepath, Knees, s,i,tag)
  % Loads all pre-computed congruity values. For each knee there are three potential visits, so three values.
  % For non-found visits, the value is set to NaN.
  if ~exist(savepath,'dir')
     error(['Could not find ',savepath])
  end
  mats = dir([savepath,'/*.mat']);
  count = length(mats);
  if ~count, pathstr, error(['No smoothness in ',savepath]), end
  kneeCount = length(Knees);
  Congruity = zeros(kneeCount,3);
  Congruity(:) = NaN;
  for m = 1:count
     file = mats(m).name;
     if strcmp(tag,'Manual')
         if isnan(str2double(file(1)))
            knee = [file(2:11),' ',file(13)]; % Becomes CCBR number and then L or R
            visit = 3; % Visit is 3 for baseline rescan
         else
            knee = [file(1:10),' ',file(12)];
            visit = 1; % Visit is 1 for basaline scan
         end
         load([savepath,'/',file]);
     % Check if the knees has already been found
        idx = strcmp(knee,Knees);
     else
        knee = [file(1:10),' ',file(30)]; % Becomes CCBR number and then L or R
        visit = str2double(file(end-4)); % Visit is 1 for baseline and 3 for baseline rescanned
        if visit==6, visit=2; end        % Visit is 6 for follow-up, we put that in the second place
        load([savepath,'/',file]);
        % Check if the knees has already been found
        idx = strcmp(knee,Knees);
        if isempty(idx)
           error(['Did not find knee: ',knee])
        end
     end
     try
         Congruity(idx,visit) = 1/(Congru(s).Cong(i).congruity.dcongc);
         %Congruity(idx,visit) = Congru(s).Cong(i).oldCong/5;
     catch
         Congruity(idx,visit) = NaN;
     end
  end

function ok = Precision(Val,tag)
  V1 = Val(:,1);
  V2 = Val(:,3);
  ok = ~isnan(V1) & ~isnan(V2);
  V1 = V1(ok);
  V2 = V2(ok);
  [c,p] = corrcoef(V1, V2);
  stds = std([V1'; V2']);
  means = mean([V1'; V2']);
  cvs = stds ./ means;
  fprintf('%s (on %d scans): RMS CV %4.1f%% Mean CV %4.1f%% Linear CC %0.2f (p %0.15f)% \n\n\n',tag,length(means),...
      100*sqrt(nanmean(cvs.^2)),100*nanmean(cvs), c(2), p(2));

  
  function voldif = longchange(nVolTF, edi, ok1)
   switch edi
       case 1
           voldif = (nVolTF(ok1,2) - nVolTF(ok1,1));
       case 2
           voldif = (nVolTF(ok1,2) - nVolTF(ok1,1))./nVolTF(ok1,1);
       case 3
           voldif = (nVolTF(ok1,2) - nVolTF(ok1,1))./Interval(ok1);
       case 4
           voldif = ((nVolTF(:,2) - nVolTF(:,1))./nVolTF(:,1)).^(1./Interval);
       case 5
           voldif = ((nVolTF(ok1,2)./nVolTF(ok1,1)).^(1./Interval(ok1))) - 1;
       otherwise
   end
   
