function test

sel = input('Choose 1 or 2 or 3:');
rows = 12;
cols = 12;
img1=zeros(rows,cols-1);
centroid(1) = round(rows/2);
centroid(2) = round(cols/2);
radius = round(rows/2) - round(rows/6);

for x=1:cols
    for y=1:rows
        r=sqrt((x-centroid(1))^2+(y-centroid(2))^2);
        r1=sqrt((x-(centroid(1)+1))^2+(y-(centroid(2)))^2);
           if r <= radius & x <= round(rows/2)
               img1(x,y)=1;
           elseif r1 <= radius & x > round(rows/2)
               img1(x,y)=2;
           end
    end
end

% Changing case
% Incongruent
img2 = img1;
start = round(rows/6);
mid = rows/2;
clo = rows - round(rows/6) + 1;
img2(start:mid,:) = img1(mid+1:clo,:);
img2(mid+1:clo,:) = img1(start:mid,:);

%Congruent
img3 = zeros(rows,cols-1);
img3(2:end,1) = 2;img3(1,1) = 1;
img3(3:end,2) = 2;img3(1:2,2) = 1;
img3(4:end,3) = 2;img3(1:3,3) = 1;
img3(5:end,4) = 2;img3(1:4,4) = 1;
img3(6:end,5) = 2;img3(1:5,5) = 1;
img3(7:end,6) = 2;img3(1:6,6) = 1;
img3(6:end,7) = 2;img3(1:5,7) = 1;
img3(5:end,8) = 2;img3(1:4,8) = 1;
img3(4:end,9) = 2;img3(1:3,9) = 1;
img3(3:end,10) = 2;img3(1:2,10) = 1;
img3(2:end,11) = 2;img3(1,11) = 1;

% Tracing the boundary

[a,b] = find(img1 == 1);
p = a(del2(b)<0);
p(end+1) = p(1);
q = unique(b);

[c,d] = find(img1 == 2);
r1 = c(del2(d)>0);
r = zeros(length(r1)+1,1);
r(1) = r1(end);r(2:end) = r1;
s = unique(d);

img=zeros(rows,cols);

if sel==1
    IMG = img1;
    
for j = 1:length(p)
    img(p(j),q(j))=1;
    img(r(j),s(j))=2;
end

elseif sel==2
    IMG = img2;
    
for j = 1:length(p)
    img(p(j)+round(rows/2)-1,q(j))=1;
    img(r(j)-round(rows/2)+1,s(j))=2;
end

else
    IMG = img3;
    
for j = 1:length(p)
    if j>1
        p(j) = r(j)-1;
            img(p(j)-4,q(j))=1;
            img(r(j),s(j))=2;
    else
            img(p(j)-4,q(j))=1;
            img(r(j),s(j))=2;
    end
end
img = rot90(img,2);
end

% To visualize

img4 = single(zeros(rows,cols,rows));
for i = 1:rows
    img4(:,:,i) = img;
end

figure
s1 = smooth3(img4 == 1);
s2 = smooth3(img4 == 2);

p1 = patch(isosurface(s1));
hold on
p2 = patch(isosurface(s2));
set(p1,'FaceColor','red','EdgeColor','None')
set(p2,'FaceColor','green','EdgeColor','None')
axis off
axis tight
view([0,70]);

% Converts into 3D

circ3 = single(zeros(rows,cols-1,rows));

for i = 1:rows
    circ3(:,:,i) = IMG;
end
IMG

% Following code Computes the congruity on for test surfaces 

qc = QuantifyCong;
classImRe = circ3;

[Tibiofemoral, Tibial, Femoral] = regularise(qc, classImRe);
th = 1;
voxelsize = [1; 1; 1];
[contactarea, Fsubs, Tsubs] = quantifycontactArea(qc, Tibiofemoral, th, voxelsize);
congruity = quantifycongruity(qc, Fsubs, Tsubs, Tibial, Femoral, 0);

