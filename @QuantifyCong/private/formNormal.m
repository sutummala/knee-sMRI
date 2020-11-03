function normal = formNormal(Grad, Gmag, subs)

% this forms unit normals for corresponding voxels in the contact region
% Grad = Gradients in x, y, z directions.
% Gmag = Gradient magnitude.
% subs = Tibial or Femoral voxel position subscripts.


normalx = (Grad.Ix)./Gmag;
normaly = (Grad.Iy)./Gmag;
normalz = (Grad.Iz)./Gmag;

normal = zeros(size(subs));

for i = 1:size(subs,1)
    
    normal(i,:) = [normalx(subs(i,1), subs(i,2), subs(i,3)), normaly(subs(i,1), subs(i,2), subs(i,3)), normalz(subs(i,1), subs(i,2), subs(i,3))];

end
