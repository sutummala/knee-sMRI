function [tibProx, femProx] = plotcontactarea(qc, Fs, Ts, TF, ca)

 % This simple function plots the contact area in the medical knee joint if the
 % inputs were as follows
 % Fs = subscripts for Femoral contact region
 % Ts = subscripts for Tibial contact region
 % TF = Tibio-Femoral binary compartment
 % ca = Contact Area in mm^2.
 
   % 1=Tibial, 2=Femoral, 3=Tibial contact region, 4=Femoral contact region  
    for i = 1:length(Ts)
        TF(Ts(i,1), Ts(i,2), Ts(i,3)) = 3;
    end
    
    for i = 1:length(Fs)
        TF(Fs(i,1), Fs(i,2), Fs(i,3)) = 4;
    end
    tibProx = TF == 3;
    femProx = TF == 4;
    
%     x = round(mean(Ts(:,1)));  y = round(mean(Ts(:,2)));  z = round(mean(Ts(:,3)));
%     % plot Tibial contact area.
%     figure
%     F = smooth3(TF == 1,'gaussian',3,1); %Tibial
%     p = patch(isosurface(F,0.1));
%     set(p,'FaceColor','cyan','EdgeColor','none','AmbientStrength',.3);
%     alpha(p, 0.1)
%     hold on
%     F = smooth3(TF == 3,'gaussian',3,1); % Replace 3 by 4 for Femoral contact area
%     p1 = patch(isosurface(F,0.1));
%     set(p1,'FaceColor',[0.6 1 0.78],'EdgeColor','none','AmbientStrength',.3);
%     alpha(p1,0.3);
%     hold on
%     F = smooth3(TF == 2, 'gaussian',3,1);% Femoral
%     p2 = patch(isosurface(F,0.1));
%     set(p2, 'FaceColor',[1 0.6 0.7], 'EdgeColor', 'none', 'AmbientStrength', .3);
%     alpha(p2, 0.2); % Transparent Femoral
%     hold on
%     quiver3(y,x,z,1,0,0,50);%X-axis
%     hold on
%     quiver3(y,x,z,0,1,0,50);%Y-axis
%     hold on
%     quiver3(y,x,z,-1,-1,0,33);%Z-axis
%     lighting none
%     axis tight off
%     view(176,6)
%     title(['\fontsize{18}Contact Area is ', num2str(round(ca)), 'mm^2']);
   
   