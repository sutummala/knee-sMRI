function gradMag = GradientMagnitude(f)

  f = single(f);
  D = 0*f; % allocation
  
  % X derivative
  D(2:end-1,:,:) = (f(3:end,:,:) - f(1:end-2,:,:))/2; % central difference
  D(1      ,:,:) =  f(2    ,:,:) - f(1      ,:,:); % forward difference
  D(end    ,:,:) =  f(end  ,:,:) - f(end-1  ,:,:); % backward difference
  gradMag = D.^2;

  % Y derivative
  D(:,2:end-1,:) = (f(:,3:end,:) - f(:,1:end-2,:))/2; % central difference
  D(:,1      ,:) =  f(:,2    ,:) - f(:,1      ,:); % forward difference
  D(:,end    ,:) =  f(:,end  ,:) - f(:,end-1  ,:); % backward difference
  gradMag = gradMag + D.^2;

  % Z derivative
  D(:,:,2:end-1) = (f(:,:,3:end) - f(:,:,1:end-2))/2; % central difference
  D(:,:,1      ) =  f(:,:,2    ) - f(:,:,1      ); % forward difference
  D(:,:,end    ) =  f(:,:,end  ) - f(:,:,end-1  ); % backward difference
  gradMag = sqrt(gradMag + D.^2);
  

