A framework for instantiating multi-linear algebra subroutines.  Designed on top of the libFLAME library (http://www.cs.utexas.edu/~flame/web/libFLAME.html).

Key features:
  * Allows for computation with symmetric tensors notably implementing sttsm operation (symmetric tensor times same matrix in all modes)
    1. Reduces storage and computational costs by employing BCSS approach to storage and algorithm design
    1. Implements variant of sttsm where matrix is multiplied in all-but-one modes
  * Windows-friendly
  * Provides a MATLAB Mex interface