function  progressionanalysis


load('overviewFull.mat','KL', 'n*','I*');
load('progressiondata.mat');
load('cong.mat');

OAKL = input('Enter the KL: ');
a = find(KL(:,1) <= OAKL);
int = Interval(a); % Time interval b/w baseline and follow-up 

%% baseline and follow-up JSW values

bas = nGap(a,1);
fol = nGap(a,2);

%% baseline and follow-up Volume values in all the compartments

vtfbas = tvolume(a,1)+volume(a,1);
vtffol = tvolume(a,2)+volume(a,2);

%% baseline and follow-up C-C (cartilage-cartilage) contact area/congruity values in all the compartments

if 0
    tfbas = ContactArea(a,1);
    tffol = ContactArea(a,2);
else
    tfbas = Congruity(a,1);
    tffol = Congruity(a,2);
end

%% JSN
% dif = ((1 + ((fol-bas)./bas)).^(1./int)) - 1;
dif = (fol - bas)./bas;
% dif = ((fol./bas).^(1./int)) - 1;

%% volume loss/ change 

% vtfdif = ((1 + ((vtffol-vtfbas)./vtfbas)).^(1./int)) - 1;
vtfdif = (vtffol - vtfbas)./vtfbas;
% vtfdif = ((vtffol./vtfbas).^(1./int)) - 1;

%% congruity loss/change 
% tfdif = ((1 + ((tffol-tfbas)./tfbas)).^(1./int)) - 1;
tfdif = (tffol- tfbas)./tfbas;
% tfdif = ((tffol./tfbas).^(1./int)) - 1;

%% NP and P for JSN 

TFNP = tfbas(dif < nanmedian(dif))';
TFP = tfbas(dif >= nanmedian(dif))';

%% NP and P for volume loss/change 

VTFNP = tfbas(vtfdif < nanmedian(vtfdif))';
VTFP = tfbas(vtfdif >= nanmedian(vtfdif))';

%% STATISTICAL ANALYSIS FOR C-C CONTACT AREA/CONGRUITY AND JSN IN ALL THE COMPARTMENTS %%

fprintf('STATISTICAL ANALYSIS FOR CONTACT AREA AND JSN IN THE MEDIAL COMPARTMENT\n');

% Tibio-Femoral

ok = ~isnan(dif(:)) & ~isnan(tfdif(:));
[R,P1] = corrcoef(dif(ok),tfdif(ok));
h = chi2gof(TFNP); h1 = chi2gof(TFP);
if ~h && ~h1
[h,pvalue] = ttest2(TFNP,TFP);
tag = 'TTEST';
else
[pvalue] = ranksum(TFNP,TFP);
tag = 'Ranksum';
end
[N] = SampleSize(TFNP,TFP);
[OR,stars] = RatioVsLow(TFNP,TFP,0.50);
auc = AUC([TFNP, TFP],[zeros(1,length(TFNP)),ones(1,length(TFP))], 0);
[dA, sddA, p] = DeLongTest(TFNP,TFP,1);
fprintf('Medial Tibio-Femoral');
fprintf('p(%s) %0.11f; SS %4.0f; OR %3.1f, %3.1f; AUC %3.3f, %0.11f; CC %3.3f, %0.11f\n', tag, pvalue, N,...
    OR(2), stars(2), auc, p, R(2), P1(2));


%% STATISTICAL ANALYSIS FOR C-C CONTACT AREA/CONGRUITY AND VOLUMESS LOSS

fprintf('STATISTICAL ANALYSIS FOR CONTACT AREA AND VOLUME LOSS IN THE MEDIAL COMPARTMENT\n');

% Tibio-Femoral 

ok1 = ~isnan(vtfdif(:)) & ~isnan(tfdif(:));
[R,P1] = corrcoef(vtfdif(ok1),tfdif(ok1));
h = chi2gof(TFNP); h1 = chi2gof(TFP);
if ~h && ~h1
[h,pvalue] = ttest2(VTFNP,VTFP);
tag = 'TTEST';
else
[pvalue] = ranksum(VTFNP,VTFP);
tag = 'Ranksum';
end
[N] = SampleSize(VTFNP,VTFP);
[OR,stars] = RatioVsLow(VTFNP,VTFP,0.50);
auc = AUC([VTFNP,VTFP],[zeros(1,length(VTFNP)),ones(1,length(VTFP))], 0);
[dA, sddA, p] = DeLongTest(VTFNP,VTFP,1);
fprintf('Medial Tibio-Femoral');
fprintf('p(%s) %0.11f; SS %4.0f; OR %3.1f, %3.1f; AUC %3.3f, %0.11f; CC %3.3f, %0.11f\n', tag, pvalue, N,...
    OR(2), stars(2), auc, p, R(2), P1(2));


