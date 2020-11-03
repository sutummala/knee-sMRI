function [binSurf, binSurfROIS, binSurfROIM, binSurfROIL] = findSurffem(im, ROIthresh)

[X, Y, Z] = size(im);
im = single(im);

if ~exist('ROItresh')
    ROIthresh = [.3,.2,.1];
end


binSurf = zeros(X,Y,Z);
% binSurfTB = []; %binSurf;
% binSurfPlus = []; %binSurf;
for k = 2:X-1
    for l = 2:Z-1
        col = squeeze(im(k,:,l));
        topSurf = max(find(col == 1));
        bottomSurf = min(find(col == 1));
        if ~isempty(bottomSurf)
            binSurf(k,bottomSurf,l) = 1;
        end
    end
end
for k = 2:X-1
    for l = 2:Y-1
        col = squeeze(im(k,l,:));
        topSurf = max(find(col == 1));
        bottomSurf = min(find(col == 1));
        if ~isempty(bottomSurf)
            binSurf(k,l,bottomSurf) = 1;
        end
    end
end
for k = 2:Y-1
    for l = 2:Z-1
        col = squeeze(im(:,k,l));
        topSurf = max(find(col == 1));
        bottomSurf = min(find(col == 1));
        if ~isempty(bottomSurf)
            binSurf(bottomSurf,k,l) = 1;
        end
    end
end

vectPos = find(binSurf(:) == 1);
[XPos, YPos, ZPos] = ndgrid(1:X, 1:Y, 1:Z);
coords = [XPos(:)'; YPos(:)'; ZPos(:)'];
coordsSurf = coords(:,vectPos);
minL = min(coordsSurf(1,:));
minW = min(coordsSurf(2,:));
minH = min(coordsSurf(3,:));
l = max(coordsSurf(1,:)) - minL;
w = max(coordsSurf(2,:)) - minW;
h = max(coordsSurf(3,:)) - minH;

% Femoral ----------- ROIS = Anterior, ROIM = Central , ROIL = Posterior

ROIS = coordsSurf(1,:)>minL+ROIthresh(3)*l & coordsSurf(1,:)<minL+(1-ROIthresh(3))*l & ...
    coordsSurf(3,:)>0 & coordsSurf(3,:)<minH+(1-ROIthresh(3))*h & coordsSurf(2,:)>(w-w/4) & coordsSurf(2,:)<=w;    % Anterior
ROIM = coordsSurf(1,:)>minL+ROIthresh(3)*l & coordsSurf(1,:)<minL+(1-ROIthresh(3))*l & ...
    coordsSurf(3,:)>0 & coordsSurf(3,:)<minH+(1-ROIthresh(3))*h & coordsSurf(2,:)>w/4 & coordsSurf(2,:)<=w-w/4;    % Medial
ROIL = coordsSurf(1,:)>minL+ROIthresh(3)*l & coordsSurf(1,:)<minL+(1-ROIthresh(3))*l & ... 
     coordsSurf(3,:)>0 & coordsSurf(3,:)<minH+(1-ROIthresh(3))*h & coordsSurf(2,:)>0 & coordsSurf(2,:)<=w/4;       % Posterior
 
 
binSurfROIS = zeros(1,X*Y*Z);
binSurfROIM = zeros(1,X*Y*Z);
binSurfROIL = zeros(1,X*Y*Z);
binSurfROIS(vectPos) = ROIS;
binSurfROIM(vectPos) = ROIM;
binSurfROIL(vectPos) = ROIL;
binSurfROIS = reshape(binSurfROIS, X, Y, Z);
binSurfROIM = reshape(binSurfROIM, X, Y, Z);
binSurfROIL = reshape(binSurfROIL, X, Y, Z);
clear vectPos

[L, nbr] = bwlabeln(binSurfROIS);
[L1, nbr1] = bwlabeln(binSurfROIM);
[L2, nbr2] = bwlabeln(binSurfROIL);
componentNbr = zeros(1,nbr); componentNbr1 = componentNbr; componentNbr2 = componentNbr;

% Find the largest connected component
% Small ROI
for w = 1:nbr
    tempComponent = find(L(:) == w);
    componentNbr(w) = length(tempComponent);
end

% Medium ROI

for w = 1:nbr1
    tempComponent = find(L1(:) == w);
    componentNbr1(w) = length(tempComponent);
end

% Large ROI

for w = 1:nbr2
    tempComponent = find(L2(:) == w);
    componentNbr2(w) = length(tempComponent);
end

% Find locations for largest connected component

cn = find(componentNbr == max(componentNbr)); cn1 = find(componentNbr1 == max(componentNbr1)); cn2 = find(componentNbr2 == max(componentNbr2));
if ~isempty(cn) && ~isempty(cn1) && ~isempty(cn2)
posLCC = find(L(:) == cn(end));
posLCC1 = find(L1(:) == cn1(end));
posLCC2 = find(L2(:) == cn2(end));
binSurfLCC = zeros(1,X*Y*Z);
binSurfLCC1 = zeros(1,X*Y*Z);
binSurfLCC2 = zeros(1,X*Y*Z);
binSurfLCC(posLCC) = 1;
binSurfLCC1(posLCC1) = 1;
binSurfLCC2(posLCC2) = 1;
binSurfROIS = reshape(binSurfLCC, X,Y,Z);
binSurfROIM = reshape(binSurfLCC1, X,Y,Z);
binSurfROIL = reshape(binSurfLCC2, X,Y,Z);
end

clear posLCC posLCC1 posLCC2 cn cn1 cn2 binSurfLCC binSurfLCC1 binSurfLCC2



