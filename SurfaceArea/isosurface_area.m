function [fout, vout, cout] = isosurface_area(varargin)
%ISOSURFACE_AREA  Computes the area of an extracted Isosurface.
%  DEVELOPED FROM THE STOCK "ISOSURFACE.M" SCRIPT
%   FV = ISOSURFACE_AREA(X,Y,Z,V,ISOVALUE) computes isosurface area for
%   data V at isosurface value ISOVALUE. Arrays (X,Y,Z) specify the points
%   at which the data V is given. The struct FV contains the faces and
%   vertices of the isosurface and can be passed directly to the PATCH
%   command.
%   
%   FV = ISOSURFACE_AREA(V,ISOVALUE) assumes [X Y Z] = meshgrid(1:N, 1:M, 1:P) 
%        where [M,N,P]=SIZE(V). 
%   
%   FV = ISOSURFACE_AREA(X,Y,Z,V) or FV = ISOSURFACE_AREA(V) selects an
%   isosurface value automatically using the histogram of the
%   data. 
%   
%   FVC = ISOSURFACE_AREA(..., COLORS) interpolates the array COLORS onto
%   the scalar field and returns the interpolated values in
%   facevertexcdata. The size of the COLORS array must be the same
%   as V. 
%   
%   FV = ISOSURFACE_AREA(..., 'noshare') does not attempt to create
%   shared vertices. This is faster, but produces a larger set of
%   vertices.  
%   
%   FV = ISOSURFACE_AREA(..., 'verbose') prints progress messages to the
%        command window as the computation progresses. 
%   
%   [F, V] = ISOSURFACE_AREA(...) or  [F, V, C] = ISOSURFACE_AREA(...)
%   returns the faces and vertices (and facevertexcdata) in
%   separate arrays instead of a struct. 
%       
%
%   See also ISONORMALS, ISOCAPS, SMOOTH3, SUBVOLUME, REDUCEVOLUME,
%            REDUCEPATCH, SHRINKFACES.

%   Copyright 1984-2004 The MathWorks, Inc.
%   Developed from ISOSURFACE.m Revision: 1.8.4.5, Date: 2005/06/21 19:37:46
%   Area Calculator developed by Daniel W. Siderius (dwsideri@gmail.com)
%   $Revision: 1.0 $ $Date: 2009/09/25$


[x y z data colors value noshare verbose] = parseargs(nargin,varargin);

if length(value)>1
  error(id('ScalarIsovalue'),'Isovalue must be a scalar.'); 
end

% Take this out when other data types are handled
data = double(data);
colors = double(colors);

[msg x y z] = xyzvcheck(x,y,z,data);  error(msg) %#ok


if ~isempty(colors) && ~isequal(size(colors), size(data))
  error(id('ColorSizeMismatch'),'COLORS array must be the same size as V.'); 
end

if isempty(value)
  value = isovalue(data);
end

[v f c] = isosurf(data, colors, value, noshare, verbose);
v = v';  %these lines transpose the data
f = f';
c = c';
if isempty(v)
  v = [];
  f = [];
  c = [];
end

if ~isempty(x) && ~isempty(v)
  sz = size(x);
  if ~(isequal(x, 1:sz(2)) && isequal(y, 1:sz(1)) && isequal(z, 1:sz(3)))
    nv(:,1) = interp3( x,  v(:,1), v(:,2), v(:,3));
    nv(:,2) = interp3( y,  v(:,1), v(:,2), v(:,3));
    nv(:,3) = interp3( z,  v(:,1), v(:,2), v(:,3));
    v = nv;
  end
end

if nargout==0
  ax = [];
  fig = get(0, 'currentfigure'); 
  if ~isempty(fig)
    ax = get(fig, 'currentaxes');
  end

  p=patch('faces', f, 'vertices', v, 'facevertexcdata', value, ...
          'facecolor', 'flat', 'edgecolor', 'none', 'userdata', value);
      
  % Register handles with m-code generator
  mcoderegister('Handles',p,'Target',p,'Name','isosurface');        
 
  if ~isempty(c)
    set(p, 'facevertexcdata', c)
  end
  if ~isempty(x)
    isonormals(x,y,z,data, p);
  else
    isonormals(data, p);
  end
  
  if isempty(ax)
    view(3); 
    camlight; lighting gouraud
  end


%Begin the added source code - DWS 2009/09/25
  [ points numaxes ] = size(f);
%  points
  area = 0.0;
  for ID=1:points
    % For each triangle tile, use Heron's Formula to calculate its area, then add it to the accumulator
    a = sqrt( (v(f(ID,1),1)-v(f(ID,2),1))^2 + (v(f(ID,1),2)-v(f(ID,2),2))^2 + (v(f(ID,1),3)-v(f(ID,2),3))^2 );
    b = sqrt( (v(f(ID,2),1)-v(f(ID,3),1))^2 + (v(f(ID,2),2)-v(f(ID,3),2))^2 + (v(f(ID,2),3)-v(f(ID,3),3))^2 );
    c = sqrt( (v(f(ID,1),1)-v(f(ID,3),1))^2 + (v(f(ID,1),2)-v(f(ID,3),2))^2 + (v(f(ID,1),3)-v(f(ID,3),3))^2 );
    s = (a+b+c)/2.0;
    myarea = sqrt( s*(s-a)*(s-b)*(s-c) );
    area = area + myarea;
  end
  area

% Plot the vertex points if desired.
%  hold on
%  for ID=1:points
%    for j=1:3
%      plot3( v(f(ID,j),1),  v(f(ID,j),2), v(f(ID,j),3), '.')
%    end
%  end

elseif nargout == 1
  [ points numaxes ] = size(f);
  area = 0.0;
  for ID=1:points
    a = sqrt( (v(f(ID,1),1)-v(f(ID,2),1))^2 + (v(f(ID,1),2)-v(f(ID,2),2))^2 + (v(f(ID,1),3)-v(f(ID,2),3))^2 );
    b = sqrt( (v(f(ID,2),1)-v(f(ID,3),1))^2 + (v(f(ID,2),2)-v(f(ID,3),2))^2 + (v(f(ID,2),3)-v(f(ID,3),3))^2 );
    c = sqrt( (v(f(ID,1),1)-v(f(ID,3),1))^2 + (v(f(ID,1),2)-v(f(ID,3),2))^2 + (v(f(ID,1),3)-v(f(ID,3),3))^2 );
    s = (a+b+c)/2.0;
    myarea = sqrt( s*(s-a)*(s-b)*(s-c) );
    area = area + myarea;
  end

  fout = area;

else
  fout = f;
  vout = v;
  if ~isempty(c)
    cout = c;
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, y, z, data, colors, value, noshare, verbose] = parseargs(nin, vargin)


x = [];
y = [];
z = [];
colors = [];
value = [];
noshare = 0;
verbose = 0;

for j = 1:2
  if nin>0
    lastarg = vargin{nin};
    if ischar(lastarg)
      if ~isempty(lastarg)
	if lower(lastarg(1))=='n' % noshare
	  noshare = 1;
	end
	if lower(lastarg(1))=='v' % verbose
	  verbose = 1;
	end
      end
      nin = nin - 1;
    end
  end
end

if nin==1                % isosurface_area(v)
  data = vargin{1};
elseif nin==2            % isosurface_area(v, isoval), isosurface_area(v, colors) 
  data = vargin{1};
  value = vargin{2};
  if isequal(size(value), size(data))
    colors = value;
    value = [];
  end
elseif nin==3            % isosurface_area(v, isoval, colors) 
  data = vargin{1};
  value = vargin{2};
  colors = vargin{3};
elseif nin==4            % isosurface_area(x,y,z,v)
  x = vargin{1};
  y = vargin{2};
  z = vargin{3};
  data = vargin{4};
elseif nin==5            % isosurface_area(x,y,z,v, isovalue), isosurface_area(x,y,z,v, colors)
  x = vargin{1};
  y = vargin{2};
  z = vargin{3};
  data = vargin{4};
  value = vargin{5};
  if isequal(size(value), size(data))
    colors = value;
    value = [];
  end
elseif nin==6            % isosurface_area(x,y,z,v, isovalue, colors)
  x = vargin{1};
  y = vargin{2};
  z = vargin{3};
  data = vargin{4};
  value = vargin{5};
  colors = vargin{6};
else
  error(id('WrongNumberOfInputs'),'Wrong number of input arguments.'); 
end

function str=id(str)
str = ['MATLAB:isosurface_area:' str];

