# Homogenization
Matlab code for 2D homogenization problems

MATLAB algorithms associated with March, Carr and Turner (2018). Code will be made available once the manuscript has gone out to review

``Homogenization`` computes the effective diffusivity of a block heterogenous domain comprised of an ``m`` by ``n`` grid by solving a homogenization problem with periodic boundary conditions. The code is applicable to problems consisting of arbitrarily-sized blocks, as long as all interfaces between blocks are aligned.

The code features two implementations of the semi-analytical solution presented by March, Carr and Turner (2018).

- ``homogenization`` implements the semi-analytical method on a domain comprised of an ``m`` by ``m`` grid of blocks, where each block is the same size.

- ``homogenization_CD`` implements the semi-analytical method on a domain comprised of an ``m`` by ``n`` grid of blocks, where blocks can be of different sizes.

``homogenization_CD``is applicable to a wider class of problems however to account for this applicability will run slower than ``homogenization``

## References

If you use ``homogenization`` and/or ``homogenization_CD``, we would appreciate that you mention it in your work by citing the following paper:

## Examples

``homogenization`` and ``homogenization_CD``are best demonstrated with some examples.

### Example A

``homogenization`` can be used to compute the effective diffusivity of any geometry that can be represented as an ``m`` by ``m`` grid of blocks

The diffusivities are represented as an array in which the diffusivity in the ``(i,j)``th block is the ``(i,j)``th entry in the diffusivity array ``D``. For example, the following geometry:

<figure><img src="https://github.com/NathanMarch/Homogenization/blob/master/Figures/Example_A.jpg" width="455"></figure>

can be represented as the matrix:

``D = [1 0.1 0.1 0.1; 1 0.1 0.1 0.1; 1 1 0.1 1; 1 1 0.1 1]``

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

``Deff = [0.3014; -0.0194; -0.0194; 0.3376]``

Computing ``Deff = homogenization(D)`` yields the effective diffusivity:

``Deff = [0.2361; 0.0000; 0.0000; 0.4224]``
