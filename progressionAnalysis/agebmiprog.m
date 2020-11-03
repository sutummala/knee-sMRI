function agebmiprog

load('overviewFull.mat');
load('progressiondata.mat');
load('ca.mat');

OAKL = input('Enter the KL: ');
a = find(KL(:,1) >= OAKL);
int = Interval(a); % Time interval b/w baseline and follow-up 

%% baseline and follow-up JSW values

bas = SmoothCurv(a,1);
fol = SmoothCurv(a,2);

%% baseline and follow-up Volume values in all the compartments

vtfbas = TSmoothCurv(a,1);
vtffol = TSmoothCurv(a,2);

%% baseline and follow-up C-C contact area values in all the compartments

tfbas = Age(a,1);
tffol = Age(a,2);

%% JSN
dif = ((fol-bas)); 

%% volume loss/ change 

vtfdif = ((vtffol-vtfbas)); 

%% contact area loss/change 
tfdif = ((tffol-tfbas)); 

%% NP and P for JSN 

TFNP = tfbas(dif < nanmedian(dif))';
TFP = tfbas(dif >= nanmedian(dif))';

%% NP and P for volume loss/change 

VTFNP = tfbas(vtfdif < nanmedian(vtfdif))';
VTFP = tfbas(vtfdif >= nanmedian(vtfdif))';

%% STATISTICAL ANALYSIS FOR C-C CONTACT AREA AND JSN IN ALL THE COMPARTMENTS %%

disp(sprintf('STATISTICAL ANALYSIS FOR CONTACT AREA AND JSN IN THE MEDIAL COMPARTMENT\n'));

% Tibio-Femoral

ok = ~isnan(dif(:)) & ~isnan(tfdif(:));
[R,P1] = corrcoef(dif(ok),tfdif(ok));
[h,pvalue] = ttest2(TFNP,TFP);
[N] = SampleSize(TFNP,TFP);
[OR,stars] = RatioVsLow(TFNP,TFP,0.50);
auc = AUC([TFNP, TFP],[zeros(1,length(TFNP)),ones(1,length(TFP))], 0);
[dA, sddA, p] = DeLongTest(TFNP,TFP,1);
disp(sprintf('Medial Tibio-Femoral'));
disp(sprintf('p %0.11f; SS %4.0f; OR %3.1f, %3.1f; AUC %3.3f, %0.11f; CC %3.3f, %0.11f\n', pvalue, N,...
    OR(2), stars(2), auc, p, R(2), P1(2)));


%% STATISTICAL ANALYSIS FOR C-C CONTACT AREA AND VOLUMESS LOSS

disp(sprintf('STATISTICAL ANALYSIS FOR CONTACT AREA AND VOLUME LOSS IN THE MEDIAL COMPARTMENT\n'));

% Tibio-Femoral 

ok1 = ~isnan(vtfdif(:)) & ~isnan(tfdif(:));
[R,P1] = corrcoef(vtfdif(ok1),tfdif(ok1));
[h,pvalue] = ttest2(VTFNP,VTFP);
[N] = SampleSize(VTFNP,VTFP);
[OR,stars] = RatioVsLow(VTFNP,VTFP,0.50);
auc = AUC([VTFNP,VTFP],[zeros(1,length(VTFNP)),ones(1,length(VTFP))], 0);
[dA, sddA, p] = DeLongTest(VTFNP,VTFP,1);
disp(sprintf('Medial Tibio-Femoral'));
disp(sprintf('p %0.11f; SS %4.0f; OR %3.1f, %3.1f; AUC %3.3f, %0.11f; CC %3.3f, %0.11f\n', pvalue, N,...
    OR(2), stars(2), auc, p, R(2), P1(2)));


