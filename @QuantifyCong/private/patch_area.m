function [area_per_layer, layer_cm] = patch_area(ptch, direction)
% [area_per_layer, layer_cm] = patch_area(ptch, direction)
%
%   ptch: matlab patch representing a thresholded volume
%   direction:  coordinate selected for per-layer area & centroids; typically
%               the long direction of the patch: its axis if it were modeled
%               as a cylinder
%
%   area_per_layer: for each "y" layer, the total of all triangle areas with
%               centroid (center of mass) at that value of y
%   layer_cm:   for each "y" layer, the x-z coordinates of the center of
%               mass of the triangle with centroid at that value of y
%
% This function provides an estimate of the surface area associated with
% an isosurface, roughly attributed to each y-increment

% Ted Ballou, 10/26/05

if nargin == 1
    direction = 1;
else 
    if nargin ~= 2
        help ptch_area
        error 'wrong argument count'
    end
end

myp = get(ptch);
nf = myp.Faces;
nv = myp.Vertices;

if size(nf,2) ~= 3
    error('ptch_area expects to get triangular faces')
end

% initialize arrays to store results per "direction" layer and identify
% non-direction 2-dimensional layer coordinates
m=max(nv,[],1);
area_per_layer = zeros(1,m(direction));
layer_cm(:,2) = zeros([m(direction),1]);
layer2d = find([1,2,3] ~= direction);

for i=1:size(nf,1)  % for each face
    v1 = nv(nf(i,1),:); % first vertex
    v2 = nv(nf(i,2),:);
    v3 = nv(nf(i,3),:);
    a=norm(v1-v2);  % length of first side
    b=norm(v1-v3);
    c=norm(v2-v3);
    p = (a + b + c)/2;  % superperimeter
%   A = sqrt(p*(p - a)*(p - b)*(p - c));   % Heron's formula for triangle area 
    
    % According to wikipedia, the following formula & computation sequence
    % can be used to prevent numerical instability in the evaluation - see
    % en.wikipedia.org/wiki/Right_triangle
    vs = sort([a,b,c]);
    a1=vs(1);b1=vs(2);c1=vs(3);
    A1=(1/4)*sqrt((a1+(b1+c1))*(c1-(a1-b1))*(c1+(a1-b1))*(a1+(b1-c1)));
    
    cm=(v1+v2+v3)/3;   % This is the centroid of the triangle, its center of mass.
    cm(direction)=round(cm(direction)); % force integer so can use as index
    
    % Accumulate weighted center per layer. Take the area of this triangle multiplied
    % by the 2 non-"direction" coordinates of the centroid, and accumulate the
    % values to be summed at the index of the "direction" coordinate
    layer_cm(cm(direction),:) = layer_cm(cm(direction),:) + A1*[cm(layer2d(1)),cm(layer2d(2))]; 
    area_per_layer(cm(direction))=area_per_layer(cm(direction))+A1; % Accumulate total area
end

for i=1:m(direction)
    if area_per_layer(i)
        % assign the 2D center of mass location for each layer
        layer_cm(i,:)=layer_cm(i,:)/area_per_layer(i);
    end
end
