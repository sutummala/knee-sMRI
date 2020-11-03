function makeFig

figure;
[X,Y] = meshgrid(-1:.1:1, -1:.1:1);
Z = 1 - exp(X.^2 + Y.^2);
surf(X,Y,Z)
hold on
% Z2 = (1-exp(-X.^2 - Y.^2));
% surf(X,Y,Z2)
[X1,Y1] = meshgrid(-0.2:.05:0.2, -0.2:.05:0.2);
Z1 = 1 - exp(X1.^2 + Y1.^2);
surf(X1,Y1,Z1)
hold on
Z2 = zeros(size(X));
surf(X, Y, Z2)