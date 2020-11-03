function [newTF, Tibial, Femoral] = regularise(qc,binaryIm, scaleK, iterations, ss)



% Crop input binary image 
cropMargin = 1;
binaryIm = crop(binaryIm, cropMargin);


TF = single(SuperSampleBinary(binaryIm, ss, 'copy'));
[Tib, Kt, Ktg, GradT, HessiT, cT, caT, ccT,cpT] = meancurvatureflow(binaryIm == 1, ss, scaleK, iterations, 'Tibial');
[Fem, Kf, Kfg, GradF, HessiF, cF, caF, ccF,cpF] = meancurvatureflow(binaryIm == 2, ss, scaleK, iterations, 'Femoral');

% Make refined or reqularised Tibio-femoral binary
tib = Tib > 0; fem = Fem > 0;
newTF = combinetibfem(TF, tib, fem);

% Tibial measures
Tibial.Tib = Tib; Tibial.Kt = Kt; Tibial.Ktg = Ktg; Tibial.GradT = GradT; Tibial.HessiT = HessiT; Tibial.cT = cT; Tibial.caT = caT;
Tibial.ccT = ccT; Tibial.cpT = cpT;

% Femoral measures
Femoral.Fem = Fem; Femoral.Kf = Kf; Femoral.Kfg = Kfg; Femoral.GradF = GradF; Femoral.HessiF = HessiF; Femoral.cF = cF; Femoral.caF = caF;
Femoral.ccF = ccF; Femoral.cpF = cpF;

function newTF = combinetibfem(TF, tib, fem)

newTF = 0*TF;
[x,y,z] = size(TF);
for i = 1:x
    for j = 1:y
        for k = 1:z
            if TF(i,j,k) == 1 & find(tib(i,j,k))
                newTF(i,j,k) = 1;
            else if TF(i,j,k) == 2 & find(fem(i,j,k))
                    newTF(i,j,k) = 2;
                else
                    continue
                end
            end
        end
    end
end



function [phi, K, Kg, Grad, Hessi, Curv, CurvA, CurvC, CurvP] = meancurvatureflow(binary, ss, scaleK, iterations, tag)

    fprintf('Regularising %s compartment\n', tag);
    
    % Parameters
    dt =0.15; 
    sc = 1; % Gaussian blurring prior to flow
    narrowB = 5;
    
    % Super sample binary segmentation before smoothing
    if ss >= 1
        binary = SuperSampleBinary(binary, ss, 'copy');
      
        % Cropping for scalen (sides divisible by 2)
        sz = size(binary);
        if mod(sz(1), 2), binary(end+1,:,:) = binary(end,:,:); end
        if mod(sz(2), 2), binary(:,end+1,:) = binary(:,end,:); end
        if mod(sz(3), 2), binary(:,:,end+1) = binary(:,:,end); end
    end

    if sc > 1
        binary = real(ifftn(scalen(fftn(binary), [sc sc sc], [0 0 0])));
        binary = binary > .5; % Volume preserved ?? (No! unwanted voxels removed)
    end

    phi = single(binary - 0.5); % Boundary is now zero-crossing


    phi = signedDist3D(phi, dt, narrowB);

    for i = 1:iterations
        fprintf(' flow %d/%d\n',i,iterations)
         % Get the curvature of phi
        K = curvatureGauss3D(phi, scaleK);
        % Update phi
        phi = phi + dt*K.* GradientMagnitude(phi); % used upwind3d before (which is slightly different)
    end
    
    phi = signedDist3D(phi, dt, narrowB);
    
       
    % Compute Curvature after reqularization
    [K, Kg, Grad, Hessi] = curvatureGauss3D(phi, scaleK);
    
    if strcmp(tag,'Tibial')
        
        [curv, curvA, curvC, curvP] = findSurftib(phi > 0);
    else
        
        [curv, curvA, curvC, curvP] = findSurffem(phi > 0);
    end
    
    Kvec = K(:);
    
    posSurf = find(curv(:).*Kvec);
    Curv = nanmean(abs(Kvec(posSurf)))*ss;
    % MeCurv = median(abs(Kvec(posSurf)))*denom;
    % VCurv = iqr(abs(Kvec(posSurf)))*denom;
    % Curv.MCurv = MCurv;
    % Curv.MeCurv = MeCurv;
    % Curv.VCurv = VCurv;

    Kv = Kvec(posSurf);

    posSurf = find(curvA(:).*Kvec);
    CurvA = nanmean(abs(Kvec(posSurf)))*ss;
    % MeCurvA = median(abs(Kvec(posSurf)))*denom;
    % VCurvA = iqr(abs(Kvec(posSurf)))*denom;
    % CurvA.MCurvA = MCurvA;
    % CurvA.MeCurvA = MeCurvA;
    % CurvA.VCurvA = VCurvA;

    posSurf = find(curvC(:).*Kvec);
    CurvC = nanmean(abs(Kvec(posSurf)))*ss;
    % MeCurvC = median(abs(Kvec(posSurf)))*denom;
    % VCurvC = iqr(abs(Kvec(posSurf)))*denom;
    % CurvC.MCurvC = MCurvC;
    % CurvC.MeCurvC = MeCurvC;
    % CurvC.VCurvC = VCurvC;

    posSurf = find(curvP(:).*Kvec);
    CurvP = nanmean(abs(Kvec(posSurf)))*ss;
    % MeCurvP = median(abs(Kvec(posSurf)))*denom;
    % VCurvP = iqr(abs(Kvec(posSurf)))*denom;
    % CurvP.MCurvP = MCurvP;
    % CurvP.MeCurvP = MeCurvP;
    % CurvP.VCurvP = VCurvP;
        
    
    
     
     
