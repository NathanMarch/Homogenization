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

