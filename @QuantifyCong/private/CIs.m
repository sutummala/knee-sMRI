function [Crd, tibc, dcongc] = CIs(Tproxsubs,Fproxsubs, Tcurv, Fcurv, Tnormal, Fnormal, TGcurv, FGcurv)


lenT = size(Tproxsubs,1);
lenF = size(Fproxsubs,1);
crd = zeros(length(lenT),1);
crd1 = zeros(length(lenF),1);
tibc = crd; femc = crd; diffc = crd; dcongc = crd; acongc = crd; dcong = crd; acong = crd; fric = crd;
normc = zeros(length(crd),3);% Vector Initialization

% tibProx to femProx    
    for i = 1:lenT
        x = Tproxsubs(i,:);% Position
        edist = 0*Fproxsubs(:,1);
        
        for j = 1:lenF
            y = Fproxsubs(j,:);% Position
            edist(j) = norm(x-y); % Euclidian Distance
        end
        
        a = find(edist == min(edist), 1, 'last'); % There may be more than one femoral voxel that is at the same distance to tibial voxel.
        voxelsubscripts = [x; Fproxsubs(a,:)];
    
              
        % Compuite Distnace between normals multiplied by curvatures
        ceqratio = max(Tcurv(i),Fcurv(a))/min(Tcurv(i), Fcurv(a));
        tibsc = (Tcurv(i)) * Tnormal(i,:); % Tibial center of curvature
        femsc = (Fcurv(a)) * Fnormal(a,:); % Femoral center of curvature
        % Normalised
        ntibsc = (Tcurv(i)) * Tnormal(i,:)/(norm(Tcurv(i) * Tnormal(i,:))); % Tibial center of curvature
        nfemsc = (Fcurv(a)) * Fnormal(a,:)/(norm(Fcurv(a) * Fnormal(a,:))); % Femoral center of curvature
          
        % Possible measures for Congruity
        
%         curvd(i) = sqrt((Tcurv(i)-TGcurv(i))^2 + (Fcurv(a)-FGcurv(a))^2);
        dcongc(i) = norm((tibsc) - (femsc));% Distnace between scaled curvature vectors
        dcong(i) = norm(Tnormal(i,:) - Fnormal(a,:));
        tibc(i) = Tcurv(i); % Tibial mean/normal curvature
        femc(i) = Fcurv(a); % Femoral mean/normal curvature
        diffc(i) = tibc(i) - femc(i); % Difference between tibial and femoral mean curvatures
        acongc(i) = acos(dot(ntibsc, nfemsc));
        fric(i) = tan(acongc(i)); % coefficient of friction? 
        acong(i) = acos(dot(Tnormal(i,:), Fnormal(a,:)));
        normc(i,:) = Tnormal(i,:) - Fnormal(a,:);
    end

tibc1 = crd1; femc1 = crd1; diffc1 = crd1; dcongc1 = crd1; acongc1 = crd1; dcong1 = crd1; acong1 = crd1; 
fric1 = crd1; normc1 = zeros(length(crd1),3); % Vector Initialization

%femProx to tibProx    
    for i = 1:lenF
        x = Fproxsubs(i,:);% Position
        edist = 0*Tproxsubs(:,1);
        
        for j = 1:lenT
            y = Tproxsubs(j,:);% Position
            edist(j) = norm(x-y); % Euclidian Distance
        end
        
        a = find(edist == min(edist), 1, 'last'); % There may be more than one femoral voxel that is at the same distance to tibial voxel.
        voxelsubscripts = [x; Tproxsubs(a,:)];
    
        %Normal vectors scaled by curvatures        
        ceqratio = max(Tcurv(a),Fcurv(i))/min(Tcurv(a), Fcurv(i));
        tibsc = (Tcurv(a)) * Tnormal(a,:); % Tibial center of curvature
        femsc = (Fcurv(i)) * Fnormal(i,:); % Femoral center of curvature
        % Normalised
        ntibsc = (Tcurv(a)) * Tnormal(a,:)/(norm(Tcurv(a) * Tnormal(a,:))); % Tibial center of curvature
        nfemsc = (Fcurv(i)) * Fnormal(i,:)/(norm(Fcurv(i) * Fnormal(i,:))); % Femoral center of curvature
          
        % Possible measures for Congruity
        
%         curvd(i) = sqrt((Tcurv(a)-TGcurv(a))^2 + (Fcurv(i)-FGcurv(i))^2);
        dcongc1(i) = norm((tibsc) - (femsc));% Distance between normal curvature vectors
        dcong1(i) = norm(Tnormal(a,:) - Fnormal(i,:));
        tibc1(i) = Tcurv(a); % Tibial mean/normal curvature
        femc1(i) = Fcurv(i); % Femoral mean/normal curvature
        diffc1(i) = tibc1(i) - femc1(i); % Difference between tibial and femoral mean curvatures
        acongc1(i) = acos(dot(ntibsc,nfemsc));
        fric1(i) = tan(acongc1(i));
        acong1(i) = acos(dot(Tnormal(a,:), Fnormal(i,:)));
        normc1(i,:) = Tnormal(a,:) - Fnormal(i,:);
     end    
% Saving the measures over Contact Area as a single structure

   nc = length(find(tibc.*femc < 0))/length(tibc); % Percentage of local pairs having the same sign
   nc1 = length(find(tibc1.*femc1 < 0))/length(tibc1); 
   Crd.dcong = mean([median(abs(dcong)),median(abs(dcong1))]);% Distance between local normal vectors
   Crd.dcongc = mean([median(abs(dcongc)),median(abs(dcongc1))]);% Distance between local curvature scaled normal vectors 
   Crd.tibc = mean([median(abs(tibc)),median(abs(tibc1))]);% Tibial proximity curvatures
   Crd.femc = mean([median(abs(femc)),median(abs(femc1))]);% Femoral proximity curvatures
   Crd.diffc = mean([median(abs(diffc)),median(abs(diffc1))]);% Distance between tibial and femoral proximity curvatures
   Crd.acong = mean([median(abs(acong)),median(abs(acong1))]);% Angle between local normal vectors
   Crd.acongc = mean([median(abs(acongc)),median(abs(acongc1))]);% Angle between local curvature scaled normal vectors
   Crd.nc = mean([nc,nc1]);
   Crd.fric = mean([var(diffc) * sum(var(normc)), var(diffc1) * sum(var(normc1))]);
   Crd.fric1 = mean([median(abs(fric)), median(abs(fric1))]);
   
   
