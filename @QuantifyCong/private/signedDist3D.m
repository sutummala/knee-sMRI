function sDist = signedDist3D(phi, dt, narrowB)
  % function sDist = signedDist(phi, dt, narrowB);
  % Originally created by: Jacob Raundahl, Sune Keller and Jenny Folkesson, 2004
  % Modified 2008,2009 by Erik Dam
  % Purpose:    Get a signed distance map of phi
  % Input parameters:   phi:        phi
  %                     dt:         Time step
  %                     narrowB:    Iterations when calculating distance map
  % ----------------------------------------------------------------------------

  % Get the distance map for phi and -phi
  dist1 = distmap(phi,dt, narrowB);
  dist2 = distmap(-phi,dt, narrowB);

  sDist = max(dist1, 0) - max(dist2,0);

  % Boundary issues
  sDist(1,:,:) = sDist(2,:,:);
  sDist(end,:,:) = sDist(end-1,:,:);
  sDist(:,1,:) = sDist(:,2,:);
  sDist(:,end,:) = sDist(:,end-1,:);
  sDist(:,:,1) = sDist(:,:,2);
  sDist(:,:,end) = sDist(:,:,end-1);
  
% ----------------------------------------------
function dist = distmap(f, dt, T)
  % dist = distmap(f, dt, T);
  % Input parameters:   f:      phi or -phi
  %                     dt:     time step
  %                     T:      Stopping criteria
  % Output parameters:  dist:   Distance map
  % ----------------------------------------------

  dist = 0*f; % allocation
  dist(f < 0) = -inf;
  dist(f > 0) = inf;

  % Calculate distance map until stopping criteria is met (narrow band)
  tp = 0;
  while tp < T*dt
     % gradMag = upwindVec3D(f);
     tp = tp + dt;
     f_new = f - dt * GradientMagnitude(f);
     dist(f > 0 & f_new < 0) = tp;
     f = f_new;
  end
 
  % Set pixels outside narrow band to plateau with min/max values
  neg = (dist == -Inf);
  pos = (dist == Inf);
  dist(neg) = 0;
  dist(pos) = 0;
  dist(neg) = min(dist(:));
  dist(pos) = max(dist(:));
