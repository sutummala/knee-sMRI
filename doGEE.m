function doGEE

% This function arrages data for GEE analysis, model the inter-knee
% Correaltions within subjects (correction for intra-knee correlations)
% 14/12-2011 by Tummala
% After arranged the data, GEEQBOX toolbox would be used to compute the regression
% parameters

load('C:\Documents and Settings\stu.SUDHAKARTUMMALA\Desktop\Congruity\Matlab\ProgressionAnalysis\overviewFull.mat');
load('C:\Documents and Settings\stu.SUDHAKARTUMMALA\Desktop\Congruity\Matlab\CACI.mat');
load('C:\Documents and Settings\stu.SUDHAKARTUMMALA\Desktop\Congruity\Matlab\dataforGEE1.mat');
addpath('C:\Documents and Settings\stu.SUDHAKARTUMMALA\Desktop\GEEQBOX');

size = length(KL(KL(:,1) == 0,:));
step = 4; % 2+2 for BL+FU

% Create subject ids (2 for baseline, 2 for follow-up for left and right
% knees)

k = 1;
subjects = ones(size*2,1);
for i = 1:step:size*2
subjects(i:i+3,1) = k;
k = k+1;
end

% Arrange data and Covariates
CA = ones(size*2,1); t = CA; side = CA; CI = CA; m = 1;
X = ones(size*2,4);

for j = 1:step:size*2
    %CA or CI values
    CA(j:j+1,1) = ContactArea(m:m+1,1);
    CA(j+2:j+3,1) = ContactArea(m:m+1,2);
    CI(j:j+1,1) = Congruity(m:m+1,1);
    CI(j+2:j+3,1) = Congruity(m:m+1,2);
    % Sex
    X(j:j+1,1) = Sex(m:m+1,1);
    X(j+2:j+3,1) = Sex(m:m+1,2);
    % Age
    X(j:j+1,2) = Age(m:m+1,1);
    X(j+2:j+3,2) = Age(m:m+1,2);
    % BMI
    X(j:j+1,3) = BMI(m:m+1,1);
    X(j+2:j+3,3) = BMI(m:m+1,2);
    %time (BL or FU)
    t(j:j+1,1) = 1;
    t(j+2:j+3,1) = 0;
    %left or right (knees)
    side(j:j+1,1) = [1,0];
    side(j+2:j+3,1) = [1,0];
    m = m+2;
end

save('dataforGEE', 'subjects', 'CA', 'CI', 'X', 't', 'side')
ok = ~isnan(CA) & ~isnan(CI) & ~isnan(X(:,2)) & ~isnan(X(:,3));
[bhat,alpha,results] = gee(subjects(ok),CI(ok),side(ok),X(ok,:));% GEE analysis using GEEQBOX toolbox

   