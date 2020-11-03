function SArea

tic
datapath = 'C:\Documents and Settings\stu\Desktop\Smoothness\Save';
segpath = [datapath,'/Segmentations'];
savepath = 'C:\Documents and Settings/stu/Desktop/Congruity/Save/surfacearea_femoral';
if ~exist(savepath,'dir'), mkdir(savepath), end
segs = dir([segpath,'/*.mat']);

for i = 1:length(segs)
    fprintf('Computing surfacearea of segmentation %3d of %3d\n',i,length(segs));
    file=segs(i).name;
    surfacearea = Surfacearea([segpath,'/',file],[savepath,'/',file]);
    CA(i) = surfacearea;
end
toc
t=toc;

function [surfacearea] = Surfacearea(segfile,savefile)

load(segfile)
clear scan classIm
Knee = crop(classImRe,1);
im = Knee == 2;
[X, Y, Z] = size(im);
binSurf = zeros(X,Y,Z);

for k = 2:X-1
    for l = 2:Z-1
        col = squeeze(im(k,:,l));
        topSurf = find(col == 1, 1, 'last' );
        bottomSurf = find(col == 1, 1 );
        if ~isempty(bottomSurf)
            binSurf(k,bottomSurf,l) = 1;
        end
    end
end
for k = 2:X-1
    for l = 2:Y-1
        col = squeeze(im(k,l,:));
        topSurf = find(col == 1, 1, 'last' );
        bottomSurf = find(col == 1, 1 );
        if ~isempty(bottomSurf)
            binSurf(k,l,bottomSurf) = 1;
        end
    end
end
for k = 2:Y-1
    for l = 2:Z-1
        col = squeeze(im(:,k,l));
        topSurf = find(col == 1, 1, 'last' );
        bottomSurf = find(col == 1, 1 );
        if ~isempty(bottomSurf)
            binSurf(bottomSurf,k,l) = 1;
        end
    end
end

%     surfacearea = length(find(binSurf))* voxelsize(1)*voxelsize(2);

    [f, v] = isosurface(binSurf, 0.5);
    
    % NormalTri contains normals at the center of the triangle faces of the isosurface.
    NormalTri=cross(v(f(:,3),:)-v(f(:,2),:),v(f(:,2),:)-v(f(:,1),:));
    NormalTri=-1*NormalTri./repmat((sqrt(NormalTri(:,1).^2+NormalTri(:,2).^2+NormalTri(:,3).^2)),1,3);
    
    % Xc Yc Zc are the coordinates of the centers of the triangle faces of
    % the isosurface 
    X = v(:,1); Y = v(:,2); Z = v(:,3);
    Xc=mean([X(f(:,1)) X(f(:,2)) X(f(:,3))],2);
    Yc=mean([Y(f(:,1)) Y(f(:,2)) Y(f(:,3))],2);
    Zc=mean([Z(f(:,1)) Z(f(:,2)) Z(f(:,3))],2);
    
    % n contains normals at the vertices
    
    n = isonormals(binSurf, v);
    
    if size(f,1) < size(v,1)
        bottom = find(NormalTri(:,3) < 0);
    else
        bottom = find(n(:,3) < 0);
    end
    a = v(f(bottom,2), :) - v(f(bottom,1), :);
    b = v(f(bottom,3), :) - v(f(bottom,1), :);
    c = cross(a,b,2);
    surfacearea = 1/2 * sum(sqrt(sum(c.^2, 2)));
    clear a b c f v n
    save(savefile,'surfacearea')
    