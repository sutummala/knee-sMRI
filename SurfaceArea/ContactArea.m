function [Knee,contactarea] = ContactArea


f = input('Please Enter The Scan Number (1 to 620):');

datapath = '../Save/Segmentations/';
mats = dir([datapath,'/*.mat']);
file = mats(f).name;
load([datapath,'/',file])

clear scan classIm
Knee = crop(classImRe,1);
Knee = SuperSampleBinary(Knee, 1, 'copy');
[x,y,z] = size(Knee);
v = zeros(size(Knee));
c = 0;
for i=1:x
    for j=1:y
        for k=1:z
            if Knee(i,j,k) == 1 && Knee(i,j,k+1) == 2
                v(i,j,:) = Knee(i,j,:);
                c = c+1;
            end
        end
    end
end
contactarea = c* voxelsize(1)*voxelsize(2);

% Tibio-Femoral
figure
TibioFemoral = smooth3(Knee,'gaussian',3,1);
p = patch(isosurface(TibioFemoral,0.1));
set(p,'FaceColor',[.6,.6,.6],'EdgeColor','none','AmbientStrength',.3);
hold on
v = smooth3(v,'gaussian',3,1);
pv = patch(isosurface(v,0.1));
set(pv,'FaceColor',[1,1,1],'EdgeColor','none','AmbientStrength',.7);


axis tight
axis off
view(2,22);
l = light;
light('Position', -get(l,'Position'))
lighting gouraud






