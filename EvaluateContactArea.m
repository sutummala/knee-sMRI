function  [ContactArea] = EvaluateContactArea
   
   if ispc
      % savepath = '../Save/Measure_all_1801_1/';
      savepath = '../Save/Measure_all_0512_1_manual/';
   else
      savepath = '/shared/OA/QuantifyCartilage/Save/contactarea/';
   end
   
   %%  [Knees, KL1, KL2] = GetKL(CCBR);
   load('C:\Documents and Settings\stu.SUDHAKARTUMMALA\Desktop\Congruity\Matlab\ProgressionAnalysis\overviewFull.mat');
   load('C:\Documents and Settings\stu.SUDHAKARTUMMALA\Desktop\Congruity\Matlab\ProgressionAnalysis\newWidth.mat');
   % normfactor = nAreaT./AreaT;
   newWidth(:,2) = newWidth(:,1); % Follow-up uses the same as Baseline
   
   [ContactArea] = LoadContactarea(savepath, Knees, 2,2,newWidth,'Manual');
   
   %% Knees used for training of the voxel classifier were excluded
    removeKnees = [];
    for i = 1:length(Knees)
        x = strcmp(Knees(i), sort(training25));
        if sum(x), removeKnees(end+1) = i; end;
    end
    
    
    ExcludedKL = find(isnan(SmoothT(:,1))); ExcludedKL = union(ExcludedKL, removeKnees'); ExcludedKL = union(ExcludedKL, 63);
    KLOA = find(KL(:,1) >= 1); KLHealthy = find(KL(:,1) == 0);
    KLHealthy = setdiff(KLHealthy, ExcludedKL); KLOA = setdiff(KLOA, ExcludedKL);
    c = cvpartition(length(KLHealthy), 'holdout', 0.5); c1 = cvpartition(length(KLOA), 'holdout', 0.5);
    Train = [KLHealthy(training(c)); KLOA(training(c1))]; Test = [KLHealthy(test(c)); KLOA(test(c1))];
    
    base = setdiff(1:318, ExcludedKL);

    ContactArea = ContactArea(base,:); KL = KL(base,:); Sex = Sex(base,:); Interval = Interval(base); nGap = nGap(base,:); Pain = Pain(base,:);
    Width = Width(base,:); men = Sex(:,1) == 0; women = Sex(:,1) == 1; nVolTF = nVolTF(base,:); Age = Age(base,:); BMI = BMI(base,:);
    
    %% Age, BMI correct-----------------------------------    
   correct = 0;
   if correct
      [CA] = correctSomething(ContactArea, KL, Age, BMI, Sex);
      separateOA(ev, CA(:,1), KL(:,1), 'ContactArea [mm]', 1, 1, 1, 1);
   end
    
    %% Longitudinal Responsiveness as SRM----------------------------
    srm = 0;
    if srm
        basfol = ~isnan(ContactArea(:,1))& ~isnan(ContactArea(:,2));
        srm = (ContactArea(basfol,2) - ContactArea(basfol,1))./Interval(basfol); %This is not SRM, don't be panic!
        lChange = 100 * srm./(ContactArea(basfol,1)./Interval(basfol));
        KL = KL(basfol,:);
        men = Sex(basfol,1) == 0; women = Sex(basfol,1) == 1;
    end
    
   %% Age, BMI correlates with Contact Area-----
   agebmi = 0;
   if agebmi
        healthy = KL(women,1) >= 0;
        age = BMI(women,1); ca = ContactArea(women,1);
        h = chi2gof(ca); % Checking whether normally distributed or not
        if h
            [c1, p1] = corr(age(healthy), ca(healthy), 'type', 'Spearman');
        else
            [c2, p2] = corr(age(healthy), ca(healthy), 'type', 'Pearson'); %rank correlation
        end
   end
       
    %% With pain-----------------
    
    pai = 0;
    if pai
        ok1 = KL(men,2) == 0; % For '0' healthy, '1' early OA and for '>1' OA
        pain = Pain(men,2); con = ContactArea(men,2);
        con1 = pain(pain ~= 0); con2 = con(pain ~= 0);
        notnan = ~isnan(con1) & ~isnan(con2);
        [c, p] = corrcoef(con1(notnan), con2(notnan));
%         [OR, stars] = RatioVsLow(con1',con2',0.50);
%         [h, p1] = ttest2(con1, con2);
    end   
    
   %% Longitudinal change
   
efficacy = 0;
   if efficacy
    okey = KL(men,1) == 2;
    edi = 2; ca = ContactArea(okey,1);
    dif = longchange(nGap(men,:), edi, okey);
    dif1 = longchange(ContactArea(men,:), edi, okey);
   
    not = ~isnan(dif) & ~isnan(dif1) & ~isinf(dif);
    h1 = chi2gof(dif); h2 = chi2gof(dif1); % Checking whether the data is normally distributed
   
    if ~h1 && ~h2
            [cor, p] = corr(dif(not), dif1(not), 'type', 'Pearson');
    else
            [cor1, p3] = corr(dif(not), dif1(not), 'type', 'Spearman');
    end
        
    ev = Evaluate;
   end
   
[ok] = Precision(ContactArea, 'Precision'); % Reproducibility on scan-rescan
   
%% Prediction of cartilage volume loss/JSN
predictOA = 0;
if predictOA
   if 1
        prog = ContactArea((dif < nanmedian(dif)),1);
        nonprog = ContactArea((dif >= nanmedian(dif)),1);
   else
        prog = dif(ca < nanmedian(ca));
        nonprog = dif(ca >= nanmedian(ca));
   end
   
   perpoints = (nanmean(prog) - nanmean(nonprog)) * 100;
   [h,p] = ttest2(prog, nonprog);
end

% Plots the Cross-sectional separation graph based on the KL index
separateOA(ev, ContactArea(:,1), KL(:,1), 'Contact Area', 1, 2, 1, 1);


function [ContactArea] = LoadContactarea(savepath, Knees, s,i,newWidth,tag)
  % Loads all pre-computed contactarea values. For each knee there are three potential visits, so three values.
  % For non-found visits, the value is set to NaN.
  if ~exist(savepath,'dir')
     error(['Could not find ',savepath])
  end
  mats = dir([savepath,'/*.mat']);
  count = length(mats);
  if ~count, pathstr, error(['No smoothness in ',savepath]), end
  kneeCount = length(Knees);
  ContactArea = zeros(kneeCount,3);
  ContactArea(:) = NaN;
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
        idx = strmatch(knee,Knees,'exact');
     else
        knee = [file(1:10),' ',file(30)]; % Becomes CCBR number and then L or R
        visit = str2double(file(end-4)); % Visit is 1 for baseline and 3 for baseline rescanned
        if visit==6, visit=2; end        % Visit is 6 for follow-up, we put that in the second place
        load([savepath,'/',file]);
        % Check if the knees has already been found
        idx = strmatch(knee,Knees,'exact');
        if isempty(idx)
            error(['Did not find knee: ',knee])
        end
     end
     try
         ContactArea(idx,visit) = (CA(s).contA(i).contactarea)/(newWidth(idx,visit).*newWidth(idx,visit));
     catch 
         ContactArea(idx,visit) = NaN;
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
  fprintf('%s (on %d scans): RMS CV %4.1f%% Mean CV %4.1f%% Linear CC %0.2f (p %0.15f)% \n\n',tag,length(means),...
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
   

