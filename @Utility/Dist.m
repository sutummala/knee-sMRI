function [FProx, TProx] = Dist


f = input('Please Enter The Scan Number (1 to 620):');

datapath = 'C:\Documents and Settings\stu\Desktop\Smoothness\Save\Segmentations';
mats = dir([datapath,'/*.mat']);
file = mats(f).name;
load([datapath,'/',file])

clear scan classIm
TF = crop(classImRe,1); 
% TF is the Tibio-Femoral binary compartment.
% Tibial voxel assigned 1 and Femoral voxel assigned 2

if size(voxelsize,2)==1
voxelsize = voxelsize';
end

wrongAnisotropy = 0;
[FProx, TProx] = CalcUnsignedDistMap(TF, voxelsize, wrongAnisotropy);
  
%   if min(voxelsize)~=max(voxelsize)
%      wrongAnisotropy = 1;
%      distMap2 = CalcUnsignedDistMap(TF, voxelsize, wrongAnisotropy);     
%      distMap = min(distMap,distMap2);
%   end
  
  
function [FProx, TProx] = CalcUnsignedDistMap(TF, voxelsize, wrongAnisotropy)
    
    size_TF = size(TF);
    
    % Perimeter of Tibial and Femoral Compartments 
    T_perim = bwperim(TF==1,conndef(ndims(TF),'maximal'));
    F_perim = bwperim(TF==2,conndef(ndims(TF),'maximal'));

    % Finds Tibial perimeter subscripts from the indices
    T_perim_idx = find(T_perim);
    T_perim_subs = cell(1,ndims(TF));
    [T_perim_subs{:}] = ind2sub(size_TF, T_perim_idx);
    T_perim_subs = [T_perim_subs{:}];

    % Finds Femoral perimeter subscripts from the indices
    F_perim_idx = find(F_perim);
    F_perim_subs = cell(1,ndims(TF));
    [F_perim_subs{:}] = ind2sub(size_TF, F_perim_idx);
    F_perim_subs = [F_perim_subs{:}];
    
    if ~wrongAnisotropy
     % HANDLE ANISOTROPY - simply by multiplying coordinates by voxelsize
     T_perim_subs = bsxfun(@times, T_perim_subs, voxelsize);
     F_perim_subs = bsxfun(@times, F_perim_subs, voxelsize);
    end

    % Constrcut a k-d tree for Tibial Perimeter subs with bucket size of 25
     Ttree = kdtree(T_perim_subs,25);

    % Using the k-d tree, find the closest Tibial voxel for each Femoral voxel
    [dist, idx] = nnsearch(Ttree,T_perim_subs',F_perim_subs');

    % This dist measures the dist between pixel centers. If we want the dist to the boundary, we have
    % to substract half a pixel size - taking anisotropy into account. So we have to project the dist
    % vector onto the voxelsize. 
    
    disp = F_perim_subs - T_perim_subs(idx,:); % list of displacement vectors
    if wrongAnisotropy
        disp = bsxfun(@times, disp, voxelsize);
    end
    dist = sqrt(sum(disp.^2,2));
    normdisp = bsxfun(@times, disp, 1./dist); % normalized displacement vectors
    halfpixel = 0.5 * bsxfun(@times, normdisp, voxelsize);
    halfpixeldist = sqrt(sum(halfpixel.^2,2));
    dist = dist - halfpixeldist;
    
    F_perim_subs = bsxfun(@rdivide, F_perim_subs, voxelsize);
    r = dist<1;
    FProx = F_perim_subs(r,:); % Femoral contact region subscripts (Femoral Proximity)

    %     distMap = zeros(size(TF),'single');
    %     distMap(F_perim_subs) = dist;
    
    Ftree = kdtree(F_perim_subs,25);
    [dist1, idx1] = nnsearch(Ftree,F_perim_subs',T_perim_subs');
    disp1 = T_perim_subs - F_perim_subs(idx1,:);
    if wrongAnisotropy
        disp1 = bsxfun(@times, disp1, voxelsize);
    end
    dist1 = sqrt(sum(disp1.^2,2));
    normdisp1 = bsxfun(@times, disp1, 1./dist1); % normalized displacement vectors
    halfpixel1 = 0.5 * bsxfun(@times, normdisp1, voxelsize);
    halfpixeldist1 = sqrt(sum(halfpixel1.^2,2));
    dist1 = dist1 - halfpixeldist1;
    
    T_perim_subs = bsxfun(@rdivide, T_perim_subs, voxelsize);
    r1 = dist1<1;
    TProx = T_perim_subs(r1,:); % Tibial contact region subscripts (Tibial Proximity)
    
    
    % Following lines used for Plotting the FProx
    for i = 1:length(FProx)
    TF(FProx(i,1),FProx(i,2),FProx(i,3)) = 3;
    end
    F = smooth3(TF == 2,'gaussian',3,1);
    p = patch(isosurface(F,0.1));
    set(p,'FaceColor',[.6,.6,.6],'EdgeColor','none','AmbientStrength',.3);
    hold on
    v = smooth3(TF == 3,'gaussian',3,1);
    pv = patch(isosurface(v,0.1));
    set(pv,'FaceColor',[1,0.4,0.6],'EdgeColor','none','AmbientStrength',.7);
    view(14,-28);
    axis off
    axis tight
    
    % Following lines used for Plotting the TProx
    for i = 1:length(TProx)
    TF(TProx(i,1),TProx(i,2),TProx(i,3)) = 4;
    end
    figure
    F = smooth3(TF == 1,'gaussian',3,1);
    p = patch(isosurface(F,0.1));
    set(p,'FaceColor',[.6,.6,.6],'EdgeColor','none','AmbientStrength',.3);
    hold on
    v = smooth3(TF == 4,'gaussian',3,1);
    pv = patch(isosurface(v,0.1));
    set(pv,'FaceColor',[1,0.4,0.6],'EdgeColor','none','AmbientStrength',.7);
    view(14,-28);
    axis off
    axis tight
