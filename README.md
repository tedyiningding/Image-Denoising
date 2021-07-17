
# Image-Denoising
A collection of MATLAB scripts and functions which achieves image denoising based on variational methods

Some codes are adapted from [Laurent Condat's work](https://lcondat.github.io/software.html).

## TV Denoising (TVdenoise.m)
TV denoising (more specifically, the ROF model [1]) finds the unique solution <img src="https://latex.codecogs.com/svg.latex?\inline&space;\large&space;\mathbf{x}^\star&space;\in&space;\mathbb{R}^{M&space;\times&space;N}" title="\large \mathbf{x}^\star \in \mathbb{R}^{M \times N}" /> which minimises

<img src="https://latex.codecogs.com/svg.latex?\large&space;\frac{1}{2}&space;\left\Vert&space;\mathbf{x}-\mathbf{y}&space;\right\Vert&space;^2&space;&plus;&space;\lambda&space;\left\Vert&space;\mathbf{x}&space;\right\Vert_{\text{TV}}" title="\large \frac{1}{2} \left\Vert \mathbf{x}-\mathbf{y} \right\Vert ^2 + \lambda \left\Vert \mathbf{x} \right\Vert_{\text{TV}}" />

where <img src="https://latex.codecogs.com/svg.latex?\inline&space;\large&space;\mathbf{y}&space;\in&space;\mathbb{R}^{M&space;\times&space;N}" title="\large \mathbf{y} \in \mathbb{R}^{M \times N}" /> is the observed image, <img src="https://latex.codecogs.com/svg.latex?\inline&space;\large&space;\lambda" title="\large \lambda" /> is a regularisation parameter that balances the two terms, and <img src="https://latex.codecogs.com/svg.latex?\inline&space;\large&space;\left\Vert&space;\cdot&space;\right\Vert_{\text{TV}}" title="\large \left\Vert \cdot \right\Vert_{\text{TV}}" /> is the total variation which is defined by

<img src="https://latex.codecogs.com/svg.latex?\large&space;\left\Vert&space;\mathbf{x}&space;\right\Vert_{\text{TV}}&space;=&space;\left\Vert&space;\mathrm{D}&space;\mathbf{x}&space;\right\Vert_{p,1}&space;=&space;\sum_{i=1,j=1}^{M,N}&space;\left|&space;\left(&space;\mathrm{D}&space;\mathbf{x}&space;\right)_{i,j}&space;\right|_p&space;=&space;\sum_{i=1,j=1}^{M,N}&space;\left(&space;\left(&space;\mathrm{D}&space;\mathbf{x}&space;\right)_{i,j,1}^p&space;&plus;&space;\left(&space;\mathrm{D}&space;\mathbf{x}&space;\right)_{i,j,2}^p&space;\right)^{1/p}" title="\large \left\Vert \mathbf{x} \right\Vert_{\text{TV}} = \left\Vert \mathrm{D} \mathbf{x} \right\Vert_{p,1} = \sum_{i=1,j=1}^{M,N} \left| \left( \mathrm{D} \mathbf{x} \right)_{i,j} \right|_p = \sum_{i=1,j=1}^{M,N} \left( \left( \mathrm{D} \mathbf{x} \right)_{i,j,1}^p + \left( \mathrm{D} \mathbf{x} \right)_{i,j,2}^p \right)^{1/p}" />

where <img src="https://latex.codecogs.com/svg.latex?\inline&space;\large&space;\mathrm{D}:&space;\mathbb{R}^{M&space;\times&space;N}&space;\rightarrow&space;\mathbb{R}^{M&space;\times&space;N&space;\times&space;2}" title="\large \mathrm{D}: \mathbb{R}^{M \times N} \rightarrow \mathbb{R}^{M \times N \times 2}" /> is the discrete gradient operator. That is, <img src="https://latex.codecogs.com/svg.latex?\inline&space;\large&space;\left\Vert&space;\cdot&space;\right\Vert_{\text{TV}}" title="\large \left\Vert \cdot \right\Vert_{\text{TV}}" /> is the <img src="https://latex.codecogs.com/svg.latex?\inline&space;\large&space;\ell_1" title="\large \ell_1" />-norm  of the <img src="https://latex.codecogs.com/svg.latex?\inline&space;\large&space;p" title="\large p" />-norm of the pixelwise image gradients [2, pp. 168]. When <img src="https://latex.codecogs.com/svg.latex?\inline&space;\large&space;p=1" title="\large p=1" /> it is called the anisotropic TV, whereas when <img src="https://latex.codecogs.com/svg.latex?\inline&space;\large&space;p=2" title="\large p=2" /> it is called the isotropic TV. The code uses the latter.

The code solves the problem using the over-relaxed Chambolle-Pock algorithm [3, Algorithm 3.1] after obtaining a saddle-point problem [2, Example 5.6]

## TGV denoising (TGVdenoise.m)
TV regularisation only promotes piecewise constant therefore the result could suffer from staircasing artefacts (as will be seen from the results shown below). 

## Results
The clear image | The noisy image
:-:|:-:
![clear](https://github.com/tedyiningding/Image-Denoising/blob/main/images/gray.png?raw=true) | ![noisy](https://github.com/tedyiningding/Image-Denoising/blob/main/images/noisy_gray.png?raw=true)
TV denoised image | TGV denoised image
![clear](https://github.com/tedyiningding/Image-Denoising/blob/main/images/TVdenoised_gray.png?raw=true) | ![noisy](https://github.com/tedyiningding/Image-Denoising/blob/main/images/TGVdenoised_gray.png?raw=true)

The table below shows some quantitative results.
|  Images | RSNR (dB) |
| - | - |
| Noisy        | 14.1727|
| TV denoised  | 25.3879|
| TGV denoised | 25.3985|

The TGV method is only marginally higher in the RSNR but there is no obvious staircasing artefacts (see the backgrou)

## References
- [1] L. I. Rudin, S. Osher, and E. Fatemi, “Nonlinear total variation based noise removal algorithms,” _Physica D_, vol. 60, no. 1–4, pp. 259–268, 1992.
- [2] A. Chambolle and T. Pock, “An introduction to continuous optimization for imaging,” _Acta Numer._, vol. 25, pp. 161–319, 2016.
- [3] L. Condat, “A primal–dual splitting method for convex optimization involving lipschitzian, proximable and linear composite terms,” _J. Optim. Theory Appl._, vol. 158, no. 2, pp. 460–479, 2013.

<!--stackedit_data:
eyJoaXN0b3J5IjpbMTU2OTE3MTk4OSwtNDAzMTE1Mjk3LC0yMD
Q3NDk0MzkxLDYwODQzNDkxMCwtMTg3MTQzMTU4NCwtODA5NTAw
MzUwLC0xNjA4NzI5NjY1LDIzNDQ3MjE0NywyMDQ5MTk1OTEyLD
E1MTI4NjE2NTUsLTIxNDczNTU4LC0xOTA4NjE0NzEyLC0yNjYy
NzQ5OTEsMTM2OTk5ODU3NSwtMTc2OTYxMTM3OSwtMTc1Nzg1OT
A5MCwtMjQ2NjE3NzgyXX0=
-->