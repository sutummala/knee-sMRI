function [Tpro, TdistMap, tp] = TProx(TF, threshold, voxelsize)


    wrongAnisotropy = 0;
    [dist, TdistMap, T_perim_subs] = CalcUnsignedDistMap(wrongAnisotropy, threshold, TF, voxelsize);
    T_perim_subs = bsxfun(@rdivide, T_perim_subs, voxelsize');
    
    if min(voxelsize) ~= max(voxelsize)
    wrongAnisotropy = 1;
    [dist1, TdistMap, T_perim_subs1] = CalcUnsignedDistMap(wrongAnisotropy, threshold, TF, voxelsize);
    dist = min(dist,dist1);
    T_perim_subs = min(T_perim_subs, T_perim_subs1);
    end
    
    r = dist < threshold; % Thrshold is set as 1 mm. Finds the tibial region that is less than 1 mm away from the femoral surface
     % Tibial contact region subscripts (Tibial Proximity)
    tp = round(T_perim_subs(r,:));
  
    for j = 1:size(tp,1)
    TF(tp(j,1),tp(j,2),tp(j,3)) = 3;
    end
    
    if isempty(tp) && threshold < 1 % this sets the area of proximity zero if no voxels found with in voxel width
        Tpro = 0;
    else
     
    
    % Tibial prox area computation from patch
    
    [f, v] = isosurface(smooth3(TF == 3, 'box'));
%     nv = isonormals(smooth3(TF == 3, 'gaussian'), v);% vertex normals
%     
    if isempty(f) | isempty(v)
        Tpro = 0;
    else
        
    nf=cross(v(f(:,3),:)-v(f(:,2),:),v(f(:,2),:)-v(f(:,1),:));
    nf= nf./repmat((sqrt(nf(:,1).^2+nf(:,2).^2+nf(:,3).^2)),1,3); %face normals
    

    ok = nf(:,3) > 0;
    a = v(f(ok,2), :) - v(f(ok,1), :);
    b = v(f(ok,3), :) - v(f(ok,1), :);
    c = cross(a,b,2);
    Tpro = 0.5 * sum(sqrt(sum(c.^2, 2))); % Area of Tibial proximity
    clear a b c f v TP T_perim_subs ok n 
    end
    end


  
function [dist, distMap, T_perim_subs] = CalcUnsignedDistMap(wrongAnisotropy, threshold, TF, voxelsize)
    
    
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
     T_perim_subs = bsxfun(@times, T_perim_subs, voxelsize');
     F_perim_subs = bsxfun(@times, F_perim_subs, voxelsize');
    end

    % Constrcut a k-d tree for Tibial Perimeter subs with bucket size of 25
    Ftree = kdtree(F_perim_subs,25);

    % Using the k-d tree, find the closest Tibial voxel for each Femoral voxel
    [dist, idx] = nnsearch(Ftree,F_perim_subs',T_perim_subs');

    % This dist measures the dist between pixel centers. If we want the dist to the boundary, we have
    % to substract half a pixel size - taking anisotropy into account. So we have to project the dist
    % vector onto the voxelsize.
    
    try
        disp = T_perim_subs - F_perim_subs(idx,:); % list of displacement vectors
    catch
        disp = T_perim_subs;
    end
    
    if wrongAnisotropy
       disp = bsxfun(@times, disp, voxelsize');
    end
    dist = sqrt(sum(disp.^2,2));
    % To compute the distance from boundary rather than from pixel center
    if threshold == 1
    normdisp = bsxfun(@times, disp, 1./dist); % normalized displacement vectors
    halfpixel = 0.5 * bsxfun(@times, normdisp, voxelsize');
    halfpixeldist = sqrt(sum(halfpixel.^2,2));
    dist = dist - halfpixeldist;
    end
    distMap = zeros(size(TF),'single');
    distMap(T_perim_idx) = dist;
    
    
    
    
    
    
    
