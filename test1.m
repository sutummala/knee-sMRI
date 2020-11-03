function test1

sel = input('Choose case: ');
le = 12;
temp = single(zeros(le,le-1));

% sine, invert sine

if sel == 1
temp(12:end,1) = 2;temp(1,1) = 1;
temp(11:end,2) = 2;temp(1:2,2) = 1;
temp(10:end,3) = 2;temp(1:3,3) = 1;
temp(9:end,4) = 2;temp(1:4,4) = 1;
temp(8:end,5) = 2;temp(1:5,5) = 1;
temp(7:end,6) = 2;temp(1:6,6) = 1;
temp(8:end,7) = 2;temp(1:5,7) = 1;
temp(9:end,8) = 2;temp(1:4,8) = 1;
temp(10:end,9) = 2;temp(1:3,9) = 1;
temp(11:end,10) = 2;temp(1:2,10) = 1;
temp(12:end,11) = 2;temp(1,11) = 1;

% sine, sine

elseif sel == 2
temp(2:end,1) = 2;temp(1,1) = 1;
temp(3:end,2) = 2;temp(1:2,2) = 1;
temp(4:end,3) = 2;temp(1:3,3) = 1;
temp(5:end,4) = 2;temp(1:4,4) = 1;
temp(6:end,5) = 2;temp(1:5,5) = 1;
temp(7:end,6) = 2;temp(1:6,6) = 1;
temp(6:end,7) = 2;temp(1:5,7) = 1;
temp(5:end,8) = 2;temp(1:4,8) = 1;
temp(4:end,9) = 2;temp(1:3,9) = 1;
temp(3:end,10) = 2;temp(1:2,10) = 1;
temp(2:end,11) = 2;temp(1,11) = 1;

% invert sine, sine

elseif sel == 3 
temp(1:6,1) = 2;temp(7:end,1) = 1;
temp(1:5,2) = 2;temp(8:end,2) = 1;
temp(1:4,3) = 2;temp(9:end,3) = 1;
temp(1:3,4) = 2;temp(10:end,4) = 1;
temp(1:2,5) = 2;temp(11:end,5) = 1;
temp(1,6) = 2;temp(12:end,6) = 1;
temp(1:2,7) = 2;temp(11:end,7) = 1;
temp(1:3,8) = 2;temp(10:end,8) = 1;
temp(1:4,9) = 2;temp(9:end,9) = 1;
temp(1:5,10) = 2;temp(8:end,10) = 1;
temp(1:6,11) = 2;temp(7:end,11) = 1;
end

sine = single(zeros(le,le-1,le));
for i = 1:12
    sine(:,:,i) = temp;
end

qc = QuantifyCong;
classImRe = sine;

[Tibiofemoral, Tibial, Femoral] = regularise(qc, classImRe);
th = 1;
voxelsize = [1; 1; 1];
[contactarea, Fsubs, Tsubs] = quantifycontactArea(qc, Tibiofemoral, th, voxelsize);
congruity = quantifycongruity(qc, Fsubs, Tsubs, Tibial, Femoral, 0);
