function distMap = UnsignedDistMap(utilO, TF, voxelsize)  
  % This is adapting the function that bwdist uses:
  % - find the perimeter of the binary volume
  % - build a kd-tree with the coordinates as features
  % - use a nearest-neighbor-search to find find nearest perimeter point
  % - compute distance from that
  % With two additions:
  % - include voxelsize
  % - correct the voxel-to-nearest-voxel distance by subtracting half-a-voxel so that the distance
  % is to the boundary between the voxels.
  
  if 0
     % EasyButNotAnisotropic - there is no way to put voxelsize into this :-(
     distMap = bwdist(TF);
     return
  end

  if size(voxelsize,2)==1
     voxelsize = voxelsize';
  end
  
  wrongAnisotropy = 0;
  distMap = CalcUnsignedDistMap(TF, voxelsize, wrongAnisotropy);
  
  if min(voxelsize)~=max(voxelsize)
     wrongAnisotropy = 1;
     distMap2 = CalcUnsignedDistMap(TF, voxelsize, wrongAnisotropy);     
     distMap = min(distMap,distMap2);
  end
  % What is this wrongAnisotropy?????
  % Normally the anisotropy needs to be done right (by using the voxelsize on the coordinates before
  % using them in nearest-neighbor below). However, with the half-voxel correction, this may not
  % always get the shortest distance to be outside. Think of a corner. In one direction, there are
  % two 1mm voxels to the nearest outside voxel center. In the other direction, there is one 2mm voxel. 
  % Distance-wise they are equally good. But the half-voxel correction is 0.5mm in one direction and
  % 1mm in the other. So the distance to the boundary is 1mm. Therefore, in some cases, the longer
  % long-voxel distance should be chosen - because it will be shorter after half-voxel correction.
  % This is done by not including voxelsize until the nearest outside voxel neighbors have been
  % found below. This will give the wrong (longer) answer is many cases, but the shorter in this
  % special case.
  % AND then using the minimal distMap values between the two ways of doing it.
  
function distMap = CalcUnsignedDistMap(TF, voxelsize, wrongAnisotropy)
  
  size_TF = size(TF);

  % Optimization: one-valued elements of TF can be closest only to
  % two-valued elements that belong to the perimeter of TF.  By reducing
  % the number of one-valued points to be searched, we can cut down on the
  % search time.
  T_perim = bwperim(TF==1,conndef(ndims(TF),'maximal'));
  F_perim = bwperim(TF==2,conndef(ndims(TF),'maximal'));

  % Find the locations of the Tibial perimeter pixels and convert that into an
  % M-by-N array perim_subs containing the locations of M points in N-space.
  T_perim_idx = find(T_perim);
  T_perim_subs = cell(1,ndims(TF));
  [T_perim_subs{:}] = ind2sub(size_TF, T_perim_idx);
  T_perim_subs = [T_perim_subs{:}];

  % Find the locations of the Femoral perimeter pixels and convert that into a
  % P-by-N array bg_subs containing the locations of P points in N-space.
  F_perim_idx = find(F_perim);
  F_perim_subs = cell(1,ndims(TF));
  [F_perim_subs{:}] = ind2sub(size_TF, F_perim_idx);
  F_perim_subs = [F_perim_subs{:}];
  
  if ~wrongAnisotropy
     % HANDLE ANISOTROPY - simply by multiplying coordinates by voxelsize
     T_perim_subs = bsxfun(@times, T_perim_subs, voxelsize);
     F_perim_subs = bsxfun(@times, F_perim_subs, voxelsize);
  end

  % From perim_subs, construct an optimized k-d tree with a bucket size of 25.
  Ttree = kdtree(T_perim_subs,25);

  % Using the k-d tree, find the closest two-valued pixel for each one-valued pixel.
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

  distMap = zeros(size(TF),'single');
  distMap(F_perim_idx) = dist;
 
  % nnsearch can also return index if we want to know the actual perimeter element


