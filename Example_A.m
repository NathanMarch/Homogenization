% This script implements Example A of the MATLAB codes implementing a
% semi-analytical solution for two-dimensional, block heterogenous media
% available at https://github.com/NathanMarch/Homogenization

%% If you use this code please city the following publication:

% Nathan G. March, Elliot J. Carr , and Ian W. Turner,
% A fast semi-analytical homogenization method for block heterogeneous
% media, Submitted to Journal of Computational Physics, 
% https://arxiv.org/abs/1812.06680


%% Example A

D = [1 0.1 0.1 0.1; 1 0.1 0.1 0.1; 1 1 0.1 1; 1 1 0.1 1]; % Diffusivity matrix
Deff = homogenization(D); % Compute effective diffusivity

% Plot geometry
[m,n] = size(D);
xcoords = 0:1/m:1;
ycoords = 1:-1/m:0;
figure;
for i = 1:m
    for j = 1:m
        xl = xcoords(j);
        xr = xcoords(j+1);
        yb = ycoords(i);
        yu = ycoords(i+1);
        x = [xl; xr; xr; xl];
        y = [yb; yb; yu; yu];
        smap = [0.4,0.4,0.4; 0.6,0.6,0.6];
        colormap(smap)
        view(2), caxis([0,1]),
        s =  patch(x,y,D(i,j));
        s.EdgeColor = 'none';        
    end
end
ax = gca;
ax.Visible = 'off';
