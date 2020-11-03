function demographicsCCBR


% Created by Sudhakar 

clear all
clc
load('C:\Documents and Settings\stu\Desktop\Congruity\Matlab\ProgressionAnalysis\overviewFull.mat');

%%
% At Baseline
men = find(Sex(:,1) == 0); women = find(Sex(:,1) == 1); % Men and Women

fprintf('Total Men and Women participated in the study were %d and %d respectively\n', length(men)/2, length(women)/2);

KL0m = find(KL(men,1) == 0); KL1m = find(KL(men,1) == 1); KL2m = find(KL(men,1) == 2); KL3m = find(KL(men,1) > 2); % KL indices for men

KL0f = find(KL(women,1) == 0); KL1f = find(KL(women,1) == 1); KL2f = find(KL(women,1) == 2); KL3f = find(KL(women,1) > 2); % KL indices for women

fprintf('Men knees\n KL0 %d\n KL1 %d\n KL2 %d\n KL3/4 %d\n', length(KL0m), length(KL1m), length(KL2m), length(KL3m));
fprintf('Women knees\n KL0 %d\n KL1 %d\n KL2 %d\n KL3/4 %d\n', length(KL0f), length(KL1f), length(KL2f), length(KL3f));

%%
computeAge(Age, KL0m, KL1m, KL2m, KL3m, KL0f, KL1f, KL2f, KL3f, 'age');
computeAge(BMI, KL0m, KL1m, KL2m, KL3m, KL0f, KL1f, KL2f, KL3f, 'BMI');

function computeAge(Age, KL0m, KL1m, KL2m, KL3m, KL0f, KL1f, KL2f, KL3f, tag)
Age0m = Age(KL0m,1); Age1m = Age(KL1m,1); Age2m = Age(KL2m,1); Age3m = Age(KL3m,1); % KL indices Age for men

menAge = ([Age0m; Age1m; Age2m; Age3m]);

fprintf('Men %s range\n %2.1f - %2.1f (%2.1f)\n', tag, min(menAge), max(menAge), mean(menAge));

Age0f = Age(KL0f,1); Age1f = Age(KL1f,1); Age2f = Age(KL2f,1); Age3f = Age(KL3f,1); % KL indices Age for women

womenAge = ([Age0f; Age1f; Age2f; Age3f]);

fprintf('Women %s range\n %2.1f - %2.1f (%2.1f)\n', tag, min(womenAge), max(womenAge), mean(womenAge));

fprintf('Men %s\n KL0 %2.1f - %2.1f (%2.1f)\n KL1 %2.1f - %2.1f (%2.1f)\n KL2 %2.1f - %2.1f (%2.1f)\n KL3/4 %2.1f - %2.1f (%2.1f)\n',...
    tag, min(Age0m), max(Age0m), mean(Age0m), min(Age1m), max(Age1m), mean(Age1m), min(Age2m), max(Age2m), mean(Age2m),...
    min(Age3m), max(Age3m), mean(Age3m));

fprintf('Women %s\n KL0 %2.1f - %2.1f (%2.1f)\n KL1 %2.1f - %2.1f (%2.1f)\n KL2 %2.1f - %2.1f (%2.1f)\n KL3/4 %2.1f - %2.1f (%2.1f)\n',...
    tag, min(Age0f), max(Age0f), mean(Age0f), min(Age1f), max(Age1f), mean(Age1f), min(Age2f), max(Age2f), mean(Age2f), min(Age3f), max(Age3f),...
    mean(Age3f));









