% This script implements the image pixelation example of the MATLAB codes
% implementing a semi-analytical solution for two-dimensional,
% block heterogenous media% available at
% https://github.com/NathanMarch/Homogenization

%% If you use this code please city the following publication:

% Nathan G. March, Elliot J. Carr, and Ian W. Turner,
% Semi-analytical solution of the homogenization boundary value problem
% for block locally-isotropic heterogeneous media, Submitted.
% https://arxiv.org/abs/1812.06680

%% Import image and convert to array
image = 'unit_cell.png'; % import image of desired geometry
I = imread(image); A = im2bw(I,graythresh(I)); % convert image to array
D = double(A); % convert array to diffusivity array
D = D(7:518,7:518); % remove entries in array corresponding to borders of image
D(~D) = 0.1; % set values of diffusivity
[m,n] = size(D); % calculate size of diffusivity array
if m~=n
    error('Diffusivity array must be square')
end

r = 16; % pixelation parameter (number of rows of blocks in pixelated image).
% Smaller values of r correspond to more pixelated geometeries

k = m/r; % number of rows of blocks to be pixelated
if round(k) ~= k
    error('Need to choose value of r such that k is an integer. r and k should be factors of m')
end

D_average = zeros(r,r); % average diffusivity for k by k grid of blocks
for l = 1:m/k
    for j = 1:n/k
        D_average(l,j) = mean(reshape(D((l-1)*k+1:l*k,(j-1)*k+1:j*k),[k^2,1]));
        if D_average(l,j) >= 0.55
            D_average(l,j) = 1;  % round average diffusivity to 0 or 1
        else
            D_average(l,j) = 0.1;
        end
    end
end

%% Plot pixelated geometry
figure;
pcolor(D_average(end:-1:1,:)); % plot geometry
smap = [0.4,0.4,0.4; 0.6,0.6,0.6];
colormap(smap)
shading flat
ax = gca;
ax.Visible = 'off';
set(gcf,'color','w');
axis square
options.Nx = 4; % choose number of asbcissas used in quadrature rule

Deff = homogenization(D_average,options);