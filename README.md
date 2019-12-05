# Homogenization
MATLAB code for 2D homogenization problems

``Homogenization`` computes the effective diffusivity of a block heterogenous domain comprised of an ``m`` by ``n`` grid by solving a homogenization problem with periodic boundary conditions. The code is applicable to problems consisting of arbitrarily-sized blocks, as long as all interfaces between blocks are aligned.

The code features two implementations of the semi-analytical solution presented by March, Carr and Turner (2018).

- ``homogenization`` implements the semi-analytical method on a domain comprised of an ``m`` by ``m`` grid of blocks, where each block is the same size.

- ``homogenization_CD`` implements the semi-analytical method on a domain comprised of an ``m`` by ``n`` grid of blocks, where blocks can be of different sizes.

``homogenization_CD``is applicable to a wider class of problems however to account for this applicability will run slower than ``homogenization``

## References

If you use ``homogenization`` and/or ``homogenization_CD``, we would appreciate that you mention it in your work by citing the following paper:

Nathan G. March, Elliot J. Carr, and Ian W. Turner,
A fast semi-analytical homogenization method for block heterogeneous media, Submitted.
https://arxiv.org/abs/1812.06680

## Examples

``homogenization`` and ``homogenization_CD``are best demonstrated with some examples. Examples A-C consist of two materials in which dark grey blocks have diffusivity ``0.1`` and light grey blocks have diffusivity ``1``. All examples can be implemented by running the relevant script file. For examples B and C you will also need to ensure that you have Example_B.mat and Example_C.mat respectively in order to load the diffusivity matrices.

### Example A

``homogenization`` can be used to compute the effective diffusivity of any geometry that can be represented as an ``m`` by ``m`` grid of blocks

The diffusivities are represented as an array in which the diffusivity in the ``(i,j)``th block is the ``(i,j)``th entry in the diffusivity array ``D``. For example, the following geometry:

<figure><img src="https://github.com/NathanMarch/Homogenization/blob/master/Figures/Example_A.jpg" width="455"></figure>

can be represented as the matrix:

``D = [1 0.1 0.1 0.1; 1 0.1 0.1 0.1; 1 1 0.1 1; 1 1 0.1 1]``.

Computing ``Deff = homogenization(D)`` yields the effective diffusivity:

``Deff = [0.2361; 0.0000; 0.0000; 0.4224]``

### Example B

While ``homogenization`` can be used with just an array of diffusivities ``D`` as an input, there are several options that can be implemented in the code. These are options affect the size of the domain, the number of abscissas used in the integrations and the number of terms taken in the summations.

The default options are:

```
x0 = 0; % Bottom left corner of domain (0,0)
xn = 1; % Top right corner of domain (1,1)
Nx = 16; % Number of abscissas used in integrations
N = 2*Nx-3; % Number of terms used in summations
  ```
  
We can change the options to:

```
x0 = 0; % Bottom left corner of domain (0,0)
xn = 10; % Top right corner of domain (10,10)
Nx = 32; % Number of abscissas used in integrations
N = 2*Nx-3; % Number of terms used in summations
options = struct('x0',x0,'xn',xn,'Nx',Nx','N',N);
Deff = homogenization(D,options);
  ```

We use the following ``10`` by ``10`` geometry:

<figure><img src="https://github.com/NathanMarch/Homogenization/blob/master/Figures/Example_B.jpg" width="455"></figure>

and compute the effective diffusivity as:

``Deff = [0.2577 -0.0159; -0.0159 0.395]``

### Example C

We can apply our homogenization scheme to larger, complex geometries. For example, we can consider the following ``50`` by ``50`` geometry:

<figure><img src="https://github.com/NathanMarch/Homogenization/blob/master/Figures/Example_C.jpg" width="455"></figure>

and compute the effective diffusivity ``Deff = homogenization(D);``

yielding:

``Deff = [0.3741 -0.0194; -0.0194 0.2532]``

### Example D

To consider more general geometries consisting of rectangular blocks, we can use ``homogenization_CD``. For example, we can consider the following  ``4`` by ``3`` diffusivity matrix:

``D = [1 0.1 0.5; 2 1.2 0.8; 1.3 2.1 0.9; 0.3 1.1 1.4];``

To use ``homogenization_CD``, we need to specify the domain of the problem as well as the coordinates of the interfaces:

```
x0 = 0; % left side of the domain
x = [1,2]; % coordinates of vertical interfaces
xn = 3; % right side of domain
y0 = 0; % top side of domain
y = [1,2,3]; % coordinates of horizontal interfaces
ym = 4; % bottom side of domain
```

The options that can be included in ``homogenization_CD`` are:
```
Nx = 16; % number of abscissas to be used in integrals along horizontal interfaces
Ny = 16; % number of abscissas to be used in integrals along vertical interfaces
N = 29; % Number of terms used in summations
options = struct('Nx',Nx','Ny',Ny,'N',N);
```

and compute the effective diffusivity ``Deff = homogenization_CD(D,x0,xn,x,y0,ym,y,options);``

yielding:

``Deff = [0.9160 -0.0455; 0.0455 0.7927]``

## Pixellation

Block geometries can be generated by importing images into MATLAB and then converting them to an array. Once the image is converted to an array of diffusivity values, our homogenization scheme can be applied to calculate the effective diffusivity.

For example, we could consider a unit cell originally presented by Carr and Turner,
Two-scale computational modelling of water flow in unsaturated soils containing
irregular-shaped inclusions., International Journal for Numerical Methods in Engineering, 98(3), pp.
157-173.

<figure><img src="https://github.com/NathanMarch/Homogenization/blob/master/Figures/unit_cell_cropped.png" width="455"></figure>

We can import this image into MATLAB and convert the image into an array of diffusivity values. The original image is 524 by 524, so we will remove some of the outer pixels to reduce it to a 512 by 512 image.


```
I = imread(image); A = im2bw(I,graythresh(I)); % convert image to array
D = double(A); % convert array to diffusivity array
D = D(7:518,7:518); % remove entries in array corresponding to borders of image
D(~D) = 0.1; % set values of diffusivity
[m,n] = size(D); % calculate size of diffusivity array
```

Now that the image has been imported and converted to an array of diffusivity values, we can choose a pixellation parameter, r, that determines the size of the geometry that will be homogenized.

``r = 16; % pixellation parameter (number of rows of blocks in pixellated image). 
% Smaller values of r correspond to more pixellated geometeries``

This means that the pixellated geometry being homogenized will be of size 16 by 16, whereas the original image was 512 by 512. Each 32 by 32 grid of pixels in the image becomes one pixel in the pixellated image. This is done by averaging the values of the 1024 pixels in the 32 by 32 grid and rounding it to the closer of 0 or 1.

``for l = 1:r
    for j = 1:r
        D_average(l,j) = mean(reshape(D((l-1)*k+1:l*k,(j-1)*k+1:j*k),[k^2,1]));
    end
end

D_average = round(D_average); % round average diffusivity to 0 or 1
D_average(~D_average) = 0.1; % set values of diffusivity``

We can then plot the pixellated geometry:

``figure;
pcolor(D_average(end:-1:1,:)); % plot geometry
colormap(smap)
shading flat``

<figure><img src="https://github.com/NathanMarch/Homogenization/blob/master/Figures/geometry_16.pdf" width="455"></figure>

The effective diffusivity can then be calculated:

``options.Nx = 4; % choose number of asbcissas used in quadrature rule
Deff = homogenization(D,options);``

yielding:

``Deff =  [0.4694 -0.0199; -0.0199 0.4445]``.

## Pixellation Loop
We can loop through different levels of pixellation by first finding the factors of ``m``:
``values = 1:ceil(m/2);
factors = [values(rem(m,values)==0) m];``

Any element of the array ``factors`` is an acceptable value for ``k`` or ``r``, so we can loop through a subset of these values:
``k_values = factors(4:6); % choose subset of factors to be used for pixellation
Nk_values = length(k_values);``

We can then loop through the values of ``k`` to generated different pixellated geometries:
<figure><img src="https://github.com/NathanMarch/Homogenization/blob/master/Figures/geometry_4.pdf" width="455"></figure>
<figure><img src="https://github.com/NathanMarch/Homogenization/blob/master/Figures/geometry_8.pdf" width="455"></figure>
<figure><img src="https://github.com/NathanMarch/Homogenization/blob/master/Figures/geometry_16.pdf" width="455"></figure>



## Installation

``Homogenization`` can be downloaded from

https://github.com/NathanMarch/Homogenization.git

After unzipping, you will need to add the directory to the MATLAB path. You can do
this via the command:
```
addpath(homogenizationroot)
```
where `homogenizationroot` is the path to the unzipped directory.

## License

See `LICENSE.md` for licensing information.



