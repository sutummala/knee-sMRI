function Cong = IIdiff(Tib, Fem, Tproxsubs, Fproxsubs)

% Created by Tummala
% Computes the Congruence of two surfaces 
% Two surfaces are congruent if they have same 1st and 2nd fundamental
% forms

lenT = size(Tproxsubs,1);
lenF = size(Fproxsubs,1);
cong = [0*min(lenT,lenF), 1];

Tu = Tib.GradT.Ix;
Tv = Tib.GradT.Iy;
Tuu = Tib.HessiT.Ixx;
Tvv = Tib.HessiT.Iyy;
Tuv = Tib.HessiT.Ixy;

Fu = Fem.GradF.Ix;
Fv = Fem.GradF.Iy;
Fuu = Fem.HessiF.Ixx;
Fvv = Fem.HessiF.Iyy;
Fuv = Fem.HessiF.Ixy;

if lenT < lenF % Instead of 'if' and writing the same thing again, write a function that does 'tibprox' ro 'femprox' and vise versa.
    
    for i = 1:lenT
        x = Tproxsubs(i,:);% Position
        edist = 0*Fproxsubs(:,1);
    
        for j = 1:lenF
            y = Fproxsubs(j,:);% Position
            edist(j) = norm(x-y); % Euclidian Distance
            clear y
        end
        a = find(edist == min(edist), 1, 'last'); % There may be more than one femoral voxel that is at the same distance to tibial voxel.
        fs = Fproxsubs(a,:);
        
        u =Tu(x(1), x(2), x(3)); TU = [1 0 u];
        v =Tv(x(1), x(2), x(3)); TV = [0 1 v];
        n = (cross(TU,TV)/norm(cross(TU,TV)))';
        uu = Tuu(x(1), x(2), x(3)); TUU = [0 0 uu];
        vv = Tvv(x(1), x(2), x(3)); TVV = [0 0 vv];
        uv = Tuv(x(1), x(2), x(3)); TUV = [0 0 uv];
        IIT = [TUU*n TUV*n; TUV*n TVV*n];
        IT = [dot(TU,TU) dot(TU,TV); dot(TU,TV) dot(TV,TV)];
        clear u v uu vv uv n
        
        u =Fu(fs(1), fs(2), fs(3)); FU = [1 0 u];
        v =Fv(fs(1), fs(2), fs(3)); FV = [0 1 v];
        n = (cross(FU,FV)/norm(cross(FU,FV)))';
        uu = Fuu(fs(1), fs(2), fs(3)); FUU = [0 0 uu];
        vv = Fvv(fs(1), fs(2), fs(3)); FVV = [0 0 vv];
        uv = Fuv(fs(1), fs(2), fs(3)); FUV = [0 0 uv];
        IIF = [FUU*n FUV*n; FUV*n FVV*n];
        IF = [dot(FU,FU) dot(FU,FV); dot(FU,FV) dot(FV,FV)];
        Idiff = IT-IF; IIdiff = IIT-IIF;
        cong(i) = sum(abs(Idiff(:))) + sum(abs(IIdiff(:)));
                
    end
else
    for i = 1:lenF
        x = Fproxsubs(i,:);% Position
        edist = 0*Tproxsubs(:,1);
    
        for j = 1:lenT
            y = Tproxsubs(j,:);% Position
            edist(j) = norm(x-y); % Euclidian Distance
            clear y
        end
        a = find(edist == min(edist), 1, 'last'); % There may be more than one femoral voxel that is at the same distance to tibial voxel.
        fs = Tproxsubs(a,:);
        
        u =Tu(fs(1), fs(2), fs(3)); TU = [1 0 u];
        v =Tv(fs(1), fs(2), fs(3)); TV = [0 1 v];
        n = (cross(TU,TV)/norm(cross(TU,TV)))';
        uu = Tuu(fs(1), fs(2), fs(3)); TUU = [0 0 uu];
        vv = Tvv(fs(1), fs(2), fs(3)); TVV = [0 0 vv];
        uv = Tuv(fs(1), fs(2), fs(3)); TUV = [0 0 uv];
        IIT = [TUU*n TUV*n; TUV*n TVV*n];
        IT = [dot(TU,TU) dot(TU,TV); dot(TU,TV) dot(TV,TV)];
        clear u v uu vv uv n
        
        u =Fu(x(1), x(2), x(3)); FU = [1 0 u];
        v =Fv(x(1), x(2), x(3)); FV = [0 1 v];
        n = (cross(FU,FV)/norm(cross(FU,FV)))';
        uu = Fuu(x(1), x(2), x(3)); FUU = [0 0 uu];
        vv = Fvv(x(1), x(2), x(3)); FVV = [0 0 vv];
        uv = Fuv(x(1), x(2), x(3)); FUV = [0 0 uv];
        IIF = [FUU*n FUV*n; FUV*n FVV*n];
        IF = [dot(FU,FU) dot(FU,FV); dot(FU,FV) dot(FV,FV)];
        Idiff = IT-IF; IIdiff = IIT-IIF;
        cong(i) = sum(abs(Idiff(:))) + sum(abs(IIdiff(:)));
    end
end

Cong.cong = 1000 * nanmean(cong);

        
        
        
        
        
        
        
    

