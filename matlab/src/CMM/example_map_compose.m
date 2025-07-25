%%% compose_maps(params.Map_stack); % TODO
Nx = 2^7;
Ny = 2^6;

x = linspace(0,1,Nx);
y = linspace(0,1,Ny);

[X,Y] = meshgrid(x,y);

%% create two maps

% map 1
Map_X = X + sin(X*2*pi).*cos(Y*(2*pi));
Map_Y = Y + 2*cos(X*2*pi);

Map1(:,:,1) = Map_X;
Map1(:,:,2) = Map_Y;

figure(2)
subplot(1,2,1)
pcolor(X,Y,Map_X); shading flat;
subplot(1,2,2)
pcolor(X,Y,Map_Y); shading flat;

%map 2;

Map_X = X + 4*cos(X*2*pi).*cos(Y*2*pi);
Map_Y = Y + 1*cos(X*2*pi);

Map2(:,:,1) = Map_X;
Map2(:,:,2) = Map_Y;

figure(3)
subplot(1,2,1)
pcolor(X,Y,Map2(:,:,1)); shading flat;
subplot(1,2,2)
pcolor(X,Y,Map2(:,:,2)); shading flat;

%% compose two maps
% Map3 = Map1 \circ Map2
Map3_exact = X + 4*cos(2*pi*X)*cos(2*pi*Y)

mint = "linear";
% map3 X-component
Map3(:,:,1) = X + interp2(X,Y,X-Map1(:,:,1),Map2(:,:,1),Map2(:,:,2),mint);
% map3 Y-component
Map3(:,:,2) = Y + interp2(X,Y,Y-Map1(:,:,2),Map2(:,:,1),Map2(:,:,2),mint);