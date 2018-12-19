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



