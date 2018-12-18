% This script implements Example D of the MATLAB codes implementing a
% semi-analytical solution for two-dimensional, block heterogenous media
% available at https://github.com/NathanMarch/Homogenization

%% If you use this code please city the following publication:

% Nathan G. March, Elliot J. Carr , and Ian W. Turner,
% A fast semi-analytical homogenization method for block heterogeneous
% media, Submitted to Journal of Computational Physics, 
% https://arxiv.org/abs/1812.06680

%% Example D

D = [1 0.1 0.5; 2 1.2 0.8; 1.3 2.1 0.9; 0.3 1.1 1.4]; % Diffusivity matrix
x0 = 0; % left side of the domain
x = [1,2]; % coordinates of vertical interfaces
xn = 3; % right side of domain
y0 = 0; % top side of domain
y = [1,2,3]; % coordinates of horizontal interfaces
ym = 4; % bottom side of domain
Nx = 16; % number of abscissas to be used in integrals along horizontal interfaces
Ny = 16; % number of abscissas to be used in integrals along vertical interfaces
N = 29; % Number of terms used in summations
options = struct('Nx',Nx','Ny',Ny,'N',N);
Deff = homogenization_CD(D,x0,xn,x,y0,ym,y,options); % Compute effective diffusivity

