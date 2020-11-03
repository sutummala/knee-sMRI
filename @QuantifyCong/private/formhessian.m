function [Hessian] = formhessian(Hessi, sub,i)

voxel = sub(i,:);

Ixx = Hessi.Ixx(voxel(1), voxel(2), voxel(3));
Ixy = Hessi.Ixy(voxel(1), voxel(2), voxel(3));
Ixz = Hessi.Ixz(voxel(1), voxel(2), voxel(3));

Iyy = Hessi.Iyy(voxel(1), voxel(2), voxel(3));
Iyz = Hessi.Iyz(voxel(1), voxel(2), voxel(3));

Izz = Hessi.Izz(voxel(1), voxel(2), voxel(3));

Hessian = [Ixx Ixy Ixz; Ixy Iyy Iyz; Ixz Iyz Izz];

