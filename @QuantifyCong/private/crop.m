function bin = crop(bin, cropMargin)
	% Crop binary
	sz = size(bin);
	xmin = sz(1)+1; xmax = 0;
	ymin = sz(2)+1; ymax = 0;
	zmin = sz(3)+1; zmax = 0;
	for z = 1:sz(3)
		slice = bin(:,:,z);
		[xs,ys] = find(slice);
		if length(xs)
			if z<zmin, zmin=z; end
			if z>zmax, zmax=z; end
			xmin = min([xs',xmin]);	xmax = max([xs',xmax]);
			ymin = min([ys',ymin]); ymax = max([ys',ymax]);
		end
	end
	xmin = max(xmin-cropMargin,1); xmax = min(xmax+cropMargin,sz(1));
	ymin = max(ymin-cropMargin,1); ymax = min(ymax+cropMargin,sz(2));
	zmin = max(zmin-cropMargin,1); zmax = min(zmax+cropMargin,sz(3));
	bin = bin(xmin:xmax, ymin:ymax, zmin:zmax);

