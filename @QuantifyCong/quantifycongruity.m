function [Crd, dcong, acong, oldCong] = quantifycongruity(qc, FS, TS, Tibial, Femoral, ss, TF, motionplane, GAA)

% FS = Femoral subscripts belongs to the Femoral proximity.
% TS = Tibial subscripts belongs to the Tibial proximity.

%If subscripts are empty i.e. no contact region in the proximity of set
%threshold

TGrad = Tibial.GradT; FGrad = Femoral.GradF;
THessian = Tibial.HessiT; FHessian = Femoral.HessiF;


% Tibial and Femoral Gradient Magnitudes
TGmag = sqrt(TGrad.Ix.^2 + TGrad.Iy.^2 + TGrad.Iz.^2);
FGmag = sqrt(FGrad.Ix.^2 + FGrad.Iy.^2 + FGrad.Iz.^2);

% TGmag = GradientMagnitude(Tibial);
% FGmag = GradientMagnitude(Femoral);

% Tibial and femoral unit surface outward normals at every point (nothing but gradients
% in this case)
Tnormal = formNormal(TGrad, TGmag, TS);
Fnormal = formNormal(FGrad, FGmag, FS);

IIform = 0;

if IIform % Compute from Second fundamental form?
    Crd = IIdiff(Tibial, Femoral, TS, FS);
else

if motionplane
    constantvec = [1 0 0]';
    tp = constantvec/norm(constantvec);

    % compute the curvature of level surfaces of tibial and femoral
    % binary compartments
    [Tcurvatures] = normalcurvature(THessian, TS, -Tnormal, TGmag, tp);
    [Fcurvatures] = normalcurvature(FHessian, FS, -Fnormal, FGmag, tp);
    
    if GAA % Computed using Ateshian et al. 1992 (for comparision)
        [oldCong] = oldmethod(THessian, FHessian, TS, FS);
    end

    % finds the corresponding femoral voxel for each tibial voxel by
    % calculating the ecludian distance. and then
    % Quantify the Center of Radii Distance as a measure of Congruity
    [Crd, dcong, acong] = CIs(TS,FS, Tcurvatures * ss, Fcurvatures * ss, -Tnormal, -Fnormal);

else
    % Mean and Gaussian Curvatures in the contact region

    Tcurvatures = zeros(length(TS),1);
    Fcurvatures = zeros(length(FS),1);

    TGcurvatures = zeros(length(TS),1);
    FGcurvatures = zeros(length(FS),1);

    TKmean = Tibial.Kt; TKGauss = Tibial.Ktg; FKmean = Femoral.Kf; FKGauss = Femoral.Kfg;

    for i = 1:length(TS)
        Tcurvatures(i) = TKmean(TS(i,1), TS(i,2), TS(i,3));
        TGcurvatures(i) = TKGauss(TS(i,1), TS(i,2), TS(i,3));
    end

    for i = 1:length(FS)
        Fcurvatures(i) = FKmean(FS(i,1), FS(i,2), FS(i,3));
        FGcurvatures(i) = FKGauss(FS(i,1), FS(i,2), FS(i,3));
    end

    [Crd, dcong, acong] = CIs(TS,FS, Tcurvatures*ss, Fcurvatures*ss, -Tnormal, -Fnormal, TGcurvatures*ss, FGcurvatures*ss); 

end
    clear TS FS Tnormal Fnormal Tibial Femoral
end

function [normK] = normalcurvature(TgH,subs,normals,gmag, tp)

    normK = zeros(length(subs), 1);

    for i = 1:size(subs,1)
        Hess = formhessian(TgH, subs, i);
%         %Shape Index (scale invarient)
%         SI = 2/pi * atan((d(3,3)+d(2,2))/d(3,3)-d(2,2));
%         %Curvedness
%         Cn = sqrt(d(3,3)^2 + d(2,2)^2);

        n = normals(i,:);%local normal
        m = null(n);% vectors in m lie in a plane perpendicular to normal plane. 
        [v,d] = eig(m'*Hess*m);
        
        v = m*v;% Principal curvature directions
        
        pcurv = - (diag(d))/gmag(subs(i,1), subs(i,2), subs(i,3));% Principal curvatures
        
        if 0 %To visualize the local motion vectors, normals over Contact Area
            quiver3(subs(i,2), subs(i,1), subs(i,3), mvec(2), mvec(1), mvec(3),10);
            hold on
            quiver3(subs(i,2), subs(i,1), subs(i,3), n(2), n(1), n(3),10);
            hold on
            quiver3(subs(i,2), subs(i,1), subs(i,3), 0, 1, 0, 10);
        end

    %Normal Curvature along local motion vector
    mvec = cross(n,tp); % Local motion vector is cross product of local normal and sagittal axis
    normK(i) = pcurv(1)*(dot(mvec, v(:,1))^2) + pcurv(2)*(1 - dot(mvec, v(:,1))^2);% Euler's Formula
    clear Hess v d
    end
% normK = normK/max(normK);

function [oldCong] = oldmethod(THessian, FHessian, TS, FS)

lenT = size(TS,1);
lenF = size(FS,1);
oldcong = zeros(length(lenT),1);
oldcong1 = zeros(length(lenF),1); % Vector Initialization

% tibProx to femProx    
    for i = 1:lenT
        x = TS(i,:);% Position
        edist = 0*FS(:,1);
        
        for j = 1:lenF
            y = FS(j,:);% Position
            edist(j) = norm(x-y); % Euclidian Distance
        end
        
        a = find(edist == min(edist), 1, 'last'); % There may be more than one femoral voxel that is at the same distance to tibial voxel.
        y = FS(a,:);
        voxelpair = [x; y];
        
        % Tibial and Femoral hessian matrices
        Thessi = formhessian(THessian, TS, i);
        Fhessi = formhessian(FHessian, FS, a);
        oldcong(i) = computeCI(Thessi, Fhessi);
    end

%femProx to tibProx    
    for i = 1:lenF
        x = FS(i,:);% Position
        edist = 0*TS(:,1);
        
        for j = 1:lenT
            y = TS(j,:);% Position
            edist(j) = norm(x-y); % Euclidian Distance
        end
        
        a = find(edist == min(edist), 1, 'last'); % There may be more than one femoral voxel that is at the same distance to tibial voxel.
        y = TS(a,:);
        voxelpair = [x; y];
        
        % Tibial and Femoral hessian matrices
        Thessi = formhessian(THessian, TS, a);
        Fhessi = formhessian(FHessian, FS, i);
        oldcong1(i) = computeCI(Thessi, Fhessi);   
    end 
    
    oldCong = mean([mean(oldcong), mean(oldcong1)]);
    
    % Function computes congruity index using method from Atheshian
function ci = computeCI(HT, HF)

    P = poly(HT); P1 = poly(HF);  % Hessian matrices
    E = roots(P); E1 = roots(P1); % Computes the eigen values
        
    % Min and max values
      ktmin = E(2);
      ktmax = E(1);
      kfmin = E1(2);
      kfmax = E1(1);
        
    % Tibial and femoral mean and difference curvatures
    Tmean = (ktmin+ktmax)/2;
    Fmean = (kfmin+kfmax)/2;
    Tdiff = ktmin-ktmax;
    Fdiff = kfmin-kfmax;
        
    % Something to choose `ALPHA`
    a = [E(1) E(2)] > [E1(1) E1(2)];
    if a(1) == a(2)
       alpha = 0;
    else
       alpha = pi/2;
    end
        
    % Final equations that leads to congruity index value
    delta = sqrt(Tdiff^2+Fdiff^2+(2*Tdiff*Fdiff*cos(2*alpha)));
    kemin = Tmean+Fmean-(delta/2);
    kemax = Tmean+Fmean+(delta/2);
    ci = sqrt((kemin^2+(kemax)^2)/2); % CI is rms value of minimum and maximum equivalent curvature

   