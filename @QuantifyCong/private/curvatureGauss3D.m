function [Km,Kg,Grad,Hessi] = curvatureGauss3D(I, scales)
% Purpose: Calculate isophote curvature of image
% Input parameters:   I:  3D Image
% Output parameters:  K:  The curvature of I

I = single(I);
ftIm = fftn(I);
clear I

scalex = scales(1); scaley = scales(2); scalez = scales(3);
% First order derivatives
Ix  = single(real(ifftn(scalen(ftIm, [scalex, scaley, scalez], [1 0 0]))));
Iy  = single(real(ifftn(scalen(ftIm, [scalex, scaley, scalez], [0 1 0]))));
Iz  = single(real(ifftn(scalen(ftIm, [scalex, scaley, scalez], [0 0 1]))));

% Second order derivatives
Ixx  = real(ifftn(scalen(ftIm, [scalex, scaley, scalez], [2 0 0])));
Iyy  = real(ifftn(scalen(ftIm, [scalex, scaley, scalez], [0 2 0])));
Izz  = real(ifftn(scalen(ftIm, [scalex, scaley, scalez], [0 0 2])));

nom = (Iyy+Izz).*(Ix.^2) + (Ixx+Izz).*(Iy.^2) + (Ixx+Iyy).*(Iz.^2);
% clear Ixx Iyy Izz

Ixy  = real(ifftn(scalen(ftIm, [scalex, scaley, scalez], [1 1 0])));
nom = nom - 2*Ix.*Iy.*Ixy;
% clear Ixy

Ixz  = real(ifftn(scalen(ftIm, [scalex, scaley, scalez], [1 0 1])));
nom = nom - 2*Ix.*Iz.*Ixz;
% clear Ixz

Iyz  = real(ifftn(scalen(ftIm, [scalex, scaley, scalez], [0 1 1])));
nom = nom - 2*Iy.*Iz.*Iyz;
Grad.Ix = Ix;
Grad.Iy = Iy;
Grad.Iz = Iz;
Hessi.Ixx = Ixx;
Hessi.Ixy = Ixy;
Hessi.Ixz = Ixz;
Hessi.Iyy = Iyy;
Hessi.Iyz = Iyz;
Hessi.Izz = Izz;
% clear Iyz

Km = double(nom./((Ix.^2 + Iy.^2 + Iz.^2).^(3/2) + eps))/2;

% % Mean curvature
% Km = ( (Iyy+Izz).*(Ix.^2) + (Ixx+Izz).*(Iy.^2) + (Ixx+Iyy).*(Iz.^2) - 2*Ix.*Iy.*Ixy - ...
%     2*Ix.*Iz.*Ixz - 2*Iy.*Iz.*Iyz )./((Ix.^2 + Iy.^2 + Iz.^2).^(3/2) + eps);

% Gaussian curvature
Kg = ( Ix.^2.*(Iyy.*Izz-Iyz.^2) + Iy.^2.*(Ixx.*Izz-Ixz.^2) + Iz.^2.*(Ixx.*Iyy-Ixy.^2) ...
    + 2*Ix.*Iy.*(Ixz.*Iyz-Ixy.*Izz) + 2*Iy.*Iz.*(Ixy.*Ixz-Iyz.*Ixx) + 2*Ix.*Iz.*(Ixy.*Iyz-Ixz.*Iyy) ) ...
    ./((Ix.^2 + Iy.^2 + Iz.^2).^2 + eps);
clear Ixx Iyy Izz Ixy Ixz Iyz

