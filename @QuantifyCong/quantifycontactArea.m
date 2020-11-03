function [contactarea, Fsubs, Tsubs] = quantifycontactArea(qc, TibioFemoral, threshold, voxelsize)


% Compute the Femoral proximity and Tibial proximity in the load area.
% contactarea is the average of the tibial and femoral proximity 


    [Fproximity, FdistMap, Fsubs] = FProx(TibioFemoral, threshold, voxelsize);
    [Tproximity, TdistMap, Tsubs] = TProx(TibioFemoral, threshold, voxelsize);
    contactarea = (Fproximity + Tproximity)/2; % contactarea in mm^2 
    clear FdistMap TdistMap