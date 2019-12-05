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

%% Calculate potential values of k and r by finding factors of m
values = 1:ceil(m/2);
factors = [values(rem(m,values)==0) m];

k_values = factors(4:6); % choose subset of factors to be used for pixelation
Nk_values = length(k_values);
Deff = zeros(2,2,Nk_values);
options.Nx = 4; % choose number of asbcissas used in quadrature rule


for i = 1:Nk_values
    k = k_values(i);
    if round(k) ~= k
        error('Need to choose value of r such that k is an integer. r and k should be factors of m')
    end
    r = m/k; % calculate pixelation parameter
    %% Plot pixelated geometry
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
    figure;
    pcolor(D_average(end:-1:1,:)); % plot geometry
    smap = [0.4,0.4,0.4; 0.6,0.6,0.6];
    colormap(smap)
    shading flat
    ax = gca;
    ax.Visible = 'off';
    set(gcf,'color','w');
    axis square
    
    Deff(:,:,i) = homogenization(D_average,options); % calculate effective diffusivity tensor
end

