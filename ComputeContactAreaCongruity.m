function ComputeContactAreaCongruity(precalc)

% maxNumCompThreads(1);
tic
if ispc
    datapath = 'C:\Documents and Settings\stu.SUDHAKARTUMMALA\Desktop\Smoothness\Save';
else
    datapath = '/nbi-fs/home/stu';
end

% Directory that contains Segmentations
segpath = [datapath,'/Segmentations'];
if ~exist(segpath,'dir'), error(['Could not find: ',segpath]), end

% Directory to save results
if ispc
    savepath = 'C:\Documents and Settings\stu.SUDHAKARTUMMALA\Desktop\Congruity\Save\Measure_all_1805_1';
else
    savepath = [datapath,'/Measure_all_2001_1'];
end
if ~exist(savepath,'dir'), mkdir(savepath), end
 
% Segmentations
 segs = dir([segpath,'/*.mat']);
if isempty(segs)
    error(['Found no Segmentations in ', segpath]);
end

% The computation loop
if ispc
    
    for s = 9:length(segs)  %KL0(s=1,2), KL1(s=5,6), KL2(s=13), KL3(s=9)
        fprintf('Computing SMOOTHNESS & CONTACTAREA & CONGRUITY for Knee%3d/%3d\n\n',s,length(segs));
        file=segs(s).name;
        DoCompute(segpath, savepath, file);
    end
    
else
    if exist('precalc','var')
         % Compute for a single scan
         precalc = str2double(precalc);
         file = segs(precalc).name;
         DoCompute(segpath, savepath, file);
    else
        for s = 1:length(segs)
            fprintf('Computing Contact Area & Fine-Scale Congruity for Knee%3d/%3d\n',s,length(segs));
            precalc = num2str(s);
            exec = ['mccexec -m 2 -j Cong',precalc,' -f ComputeContactAreaCongruity ',precalc];
            [notOK,out] = system(exec);
            if notOK
                disp(out)
                error(['Error submitting: ',exec])
            end
        end
    end
end
toc    


function DoCompute(segpath, savepath, file)
        
        qc = QuantifyCong;
        
        if exist([savepath,'/',file], 'file')
            disp(['Already Computed\n', [savepath,'/',file]])
            return
        end
        
        load([segpath,'/',file]);
        
        Congru = struct([]); TSmooth = Congru; FSmooth = Congru; TibS = Congru; FemS = Congru; Cong = Congru; contA = Congru; CA = Congru;
        
        superSample = 1;
        % Some plot/threshold settings
        threshold = 0.71; % This threshold is used to compute the tibial and femoral proximity that are voxel width apart
        areamap = 0; % Plots the contactarea if 1
        congmap = 1; % Plots the congruity if 1

               
        for i = 1:4 %Parameter
            scaleK = (1.5^(i+1) * voxelsize);
            fprintf('---------------------Curvature scale %0.3f------------------------\n', scaleK(1));
            for j = 1:4 %Parameter
                iterations = 2^(j);
                fprintf('---------%d Iterations------------\n', iterations);
                % -------------Regularise and Compute Smoothness-----------
                fprintf('Regularising the binary cartilages for knee %s\n', file(1:end-4));
                try % Automatic
                    [Tibiofemoral, Tibial, Femoral] = regularise(qc, classImRe, scaleK, iterations, superSample);
                catch % Manual 
                    [Tibiofemoral, Tibial, Femoral] = regularise(qc, binaryIm, scaleK, iterations, superSample);
                end
                %Save Tibial and Tibial ROIS curvature/smoothness
                curv = Tibial.cT; curva = Tibial.caT; curvc = Tibial.ccT; curvp = Tibial.cpT;
                TSmooth(j).curv = curv; TSmooth(j).curva = curva; TSmooth(j).curvc = curvc; TSmooth(j).curvp = curvp; 
                clear curv curva curvc curvp
                
                % Sace Tibial and Tibial anterior, central, posterior
                % curvature/smoothness
                curv = Femoral.cF; curva = Femoral.caF; curvc = Femoral.ccF; curvp = Femoral.cpF;
                FSmooth(j).curv = curv; FSmooth(j).curva = curva; FSmooth(j).curvc = curvc; FSmooth(j).curvp = curvp; 
                
                % ---------------Contact Area Computation-------------------
                fprintf('\nComputing CONTACTAREA (CA)');
                [contactarea, Fsubs, Tsubs] = quantifycontactArea(qc, Tibiofemoral, threshold, voxelsize);
                contactarea = contactarea/(superSample^2);
                if isempty(Tsubs) | isempty(Fsubs), return, end
                fprintf(' ---> %3.0fmm^2\n', contactarea);
                
                % Plot Contact Area
                if areamap
                   [tibProx, femProx] = plotcontactarea(qc, Fsubs, Tsubs, Tibiofemoral, contactarea); % Plots contactarea
                end
                
                contA(j).contactarea = contactarea;
                
                % --------------Congruity Computation----------------
                fprintf('\nComputing CONGRUITY INDEX (CI)');
                try
                    [congruity, dcong, dcongc, oldCong] = quantifycongruity(qc, Fsubs, Tsubs, Tibial, Femoral, superSample, Tibiofemoral, 1, 1);
                catch
                    [congruity] = quantifycongruity(qc, Fsubs, Tsubs, Tibial, Femoral, superSample, Tibiofemoral, 1, 1); 
                end
                fprintf(' ---> %1.2f, %1.2f\n\n', congruity.dcong, congruity.dcongc);
                %fprintf('Atheshian CI ---> %2.1fm^-1\n\n', 1000 * oldCong);
                
                % Plot Congruity/Incongruity map
                if congmap
                    congruitymap(qc, (abs(dcong)), Tsubs, Fsubs, Tibiofemoral,0);
                end
                Cong(j).congruity = congruity; % Several metrics, some unitless
                Cong(j).oldCong = 1000 * oldCong; % Its in inverse meters
            end
             clear Tibiofemoral Tibial Femoral
%              
%              % Save all combinations measure as a structure 
             TibS(i).TSmooth = TSmooth; % Tibial Smoothness/Curvature
             FemS(i).FSmooth = FSmooth; % Femoral Smoothness/Curvature
             CA(i).contA = contA; % ContactArea
             Congru(i).Cong = Cong; % Congruity
        end
    % 
        save([savepath,'/',file], 'TibS', 'FemS', 'CA', 'Congru')
%         save([savepath,'/',file], 'tibProx', 'femProx','voxelsize')
          %Ssave([savepath,'/',file], 'tibial', 'femoral', 'voxelsize')
   
 
