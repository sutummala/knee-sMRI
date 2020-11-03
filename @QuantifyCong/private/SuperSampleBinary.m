function super = SuperSampleBinary(seg, factor, method)
	% Methods:
	% - Do nothing
	% - Copy voxel to factor x factor x factors voxels
	% - Linear interpolation
	% - Cubic interpolation
	% - Spline interpolation
	% - Morphology: Circular Dilation+Erosion+erosion+dilation
	% - Gaussian blurring
	% Morph, Gaussian assumes that copy is done first
	if isempty(strmatch(method,{'copy','linear','cubic','spline','morph','gauss'}))
		super = seg;
		warning(['Unknown method: ',method])
	end

	if ~isempty(strmatch(method,{'copy','morph','gauss'}))
		% Copy voxel to factor x factor x factors voxels
		super = copyVoxel(seg, factor);
		% Continue with morphology or Gaussian if desired
		if ~isempty(strmatch(method,{'morph','gauss'}))
			error('I did not make morph or Gaussian yet!')
			% this is because it is different from in expandSlicesZ because here it actually needs to
			% be 3D !
			super = super >= 0.5;
		end
	else
		error('interp3 is just too bloody time and memory-consuming - so dont use it (at least not for factor 5)')
		% Coordinates for interp3 method:
		sz = size(seg);
		xsz = sz(1); ysz = sz(2); zsz = sz(2);
		outside = (factor-1)/(2*factor);
		xs = linspace(1-outside, xsz+outside, xsz*factor);
		ys = linspace(1-outside, ysz+outside, ysz*factor);
		zs = linspace(1-outside, zsz+outside, zsz*factor);
		[xi,yi,zi] = meshgrid(ys,xs,zs);
		outsideValue = 0;
		super = interp3(seg, xi, yi, zi, interpmethod, outsideValue);
		super = super >= 0.5;
	end
	
function super = copyVoxel(seg, factor)
	dim = size(seg);
	super = zeros(factor * dim);
	out = (factor-1)/2; % assumes that factor is odd!
	for x = 1:dim(1)
		sx = (x-1)*factor;
		for y = 1:dim(2)
			sy = (y-1)*factor;
			for z = 1:dim(3)
				sz = (z-1)*factor;
				val = seg(x,y,z);
				% Now insert in in cube at that position
				for ix = 1:factor
					for iy = 1:factor
						for iz = 1:factor
							super(sx+ix,sy+iy,sz+iz) = val;
						end
					end
				end
			end
		end
	end
