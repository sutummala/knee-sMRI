function comboMarker

load('C:\Documents and Settings\stu\Desktop\Congruity\Matlab\ProgressionAnalysis\overviewFull.mat');

% Lets develop a combintion marker for tibial and femoral volume using LDA
a = ~isnan(nVolT(:,1));
group = {'Tibial Volume', 'Femoral Volume'};
fprintf('Number of samples in each class %d, Number of classes %d \n', length(find(a)), length(group));

% Start analysis

volt = nVolT(a,1);
class1 = [(1:size(volt,1))', volt];
meant = mean(class1);

volf = nVolF(a,1);
class2 = [(1:size(volf,1))', volf];
meanf = mean(class2);

fprintf('Class means calculated\n');
% LDA is to fina a vector that maximize the between class variance to within calss variance