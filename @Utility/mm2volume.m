function points = mm2volume(utilO, points, voxelsize)
  % Rescale 
  points = bsxfun(@times, points, 1./voxelsize);
