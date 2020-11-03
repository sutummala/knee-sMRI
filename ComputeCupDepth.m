function CD = ComputeCupDepth

tic
datapath = 'C:\Documents and Settings\stu.SUDHAKARTUMMALA\Desktop\Smoothness\Save';
segpath = [datapath,'/Segmentations'];
savepath = 'C:\Documents and Settings\stu\Desktop\Congruity\Save\cupdepth_2910';
if ~exist(savepath,'dir'), mkdir(savepath), end
segs = dir([segpath,'/*.mat']);
CD = zeros(length(segs),1);

for i = 1:length(segs)
    fprintf('Computing CupDepth of segmentation %3d of %3d\n',i,length(segs));
    file=segs(i).name;
    cupdepth = CupDepth([segpath,'/',file],[savepath,'/',file]);
    fprintf('CupDepth of segmentation %3d is %1.1fmm\n', i, cupdepth);
    CD(i) = cupdepth;
end
toc
t=toc;

function [cupdepth] = CupDepth(segfile,savefile)

load(segfile)
clear scan classIm
scaleK = 1.5^2*(voxelsize); iterations = 4; superSample = 1;
qc = QuantifyCong;
[Tibiofemoral, Tibial, Femoral] = regularise(qc, classImRe, scaleK, iterations, superSample);
Tibial = crop(Tibiofemoral == 1, 1);
[x,y,z] = size(Tibial);
if mod(x,2)
    x = x+1;
end
left = zeros(z,1);
right = zeros(z,1);
top = 0;
% Find peak
for k = 1:z
    if find(Tibial(1:x/2,:,k))
        left(k) = k;
    end
    if find(Tibial(x/2+1:x-1,:,k))
        right(k) = k;
    end
end
top = min(max(left),max(right));

% Find bottom

x1 = round(0.1*x);
x2 = round(0.9*x);
y1 = round(0.1*y);
y2 = round(0.9*y);

% Intialization of sagittal and coronal axis
sag = zeros(x,1);
cor = zeros(y,1);

 for j = 1:y
     for i = 40
     [r,c] = find(Tibial(i,j,:));
     if ~isempty(c)
         sag(i) = max(c);
     end
     clear r c 
     end
     if ~isinf(1./sag)
         cor(j) = min(sag(find(sag)));
     end
 end
 bottom = min(cor(find(cor)));

cupdepth = (top-bottom)*voxelsize(2);

save(savefile,'cupdepth')

