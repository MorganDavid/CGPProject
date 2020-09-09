# Cartesian Genetic Programming (CGP) project. 

This project uses CGP and the Fourier transform to run machine learning on periodic functions with the hope of future expansion to allow predictions on real-world periodic phenomena such as energy consumption predictions.

# Code overview
Code starts in CGPProject/src. The MyFourierClass object abstracts away the work of applying the Fourier transform and decomposing into constintuent waves, simply load data, then run execute_fourier_synthesis and the decomposed waves are accesible using the class getters. 
We make use of the [FFTW library](http://www.fftw.org/) for computing the Fourier transform and the [CGPLibrary library](http://www.cgplibrary.co.uk/files2/About-txt.html) for CGP operations. 

Significant changes are made to the CGPLibrary to allow for the intergration of our system, they can be seen in cgp.c and cgp.h in [this gitdiff](https://github.com/morgandavid/cgpproject/compare/b439d5804dbdb19cf9ca698fe6f7f6197427ca2c...0a497207194d6e3a843034cf156db9d094e058cc).

MATLAB is used for experimenation and presenting results.
