
# Image-Denoising
A collection of MATLAB scripts and functions which achieves image denoising based on variational methods

Some codes are adapted from [Laurent Condat's work](https://lcondat.github.io/software.html).

## TV Denoising (TVdenoise.m)
TV denoising finds the unique solution <img src="https://latex.codecogs.com/svg.latex?\inline&space;\large&space;\mathbf{x}^\star&space;\in&space;\mathbb{R}^{M&space;\times&space;N}" title="\large \mathbf{x}^\star \in \mathbb{R}^{M \times N}" /> which minimises

<img src="https://latex.codecogs.com/svg.latex?\large&space;\frac{1}{2}&space;\left\Vert&space;\mathbf{x}-\mathbf{y}&space;\right\Vert&space;^2&space;&plus;&space;\lambda&space;\left\Vert&space;\mathbf{x}&space;\right\Vert_{\text{TV}}" title="\large \frac{1}{2} \left\Vert \mathbf{x}-\mathbf{y} \right\Vert ^2 + \lambda \left\Vert \mathbf{x} \right\Vert_{\text{TV}}" />

where <img src="https://latex.codecogs.com/svg.latex?\inline&space;\large&space;\mathbf{y}&space;\in&space;\mathbb{R}^{M&space;\times&space;N}" title="\large \mathbf{y} \in \mathbb{R}^{M \times N}" /> is the observed image, <img src="https://latex.codecogs.com/svg.latex?\inline&space;\large&space;\lambda" title="\large \lambda" /> is a regularisation parameter that balances the two terms, and <img src="https://latex.codecogs.com/svg.latex?\inline&space;\large&space;\left\Vert&space;\cdot&space;\right\Vert_{\text{TV}}" title="\large \left\Vert \cdot \right\Vert_{\text{TV}}" /> is the total variation which is defined by

<img src="https://latex.codecogs.com/svg.latex?\large&space;\left\Vert&space;\mathbf{x}&space;\right\Vert_{\text{TV}}&space;=&space;\left\Vert&space;\mathrm{D}&space;\mathbf{x}&space;\right\Vert_{p,1}&space;=&space;\sum_{i=1,j=1}^{M,N}&space;\left|&space;\left(&space;\mathrm{D}&space;\mathbf{x}&space;\right)_{i,j}&space;\right|_p&space;=&space;\sum_{i=1,j=1}^{M,N}&space;\left(&space;\left(&space;\mathrm{D}&space;\mathbf{x}&space;\right)_{i,j,1}^p&space;&plus;&space;\left(&space;\mathrm{D}&space;\mathbf{x}&space;\right)_{i,j,2}^p&space;\right)^{1/p}" title="\large \left\Vert \mathbf{x} \right\Vert_{\text{TV}} = \left\Vert \mathrm{D} \mathbf{x} \right\Vert_{p,1} = \sum_{i=1,j=1}^{M,N} \left| \left( \mathrm{D} \mathbf{x} \right)_{i,j} \right|_p = \sum_{i=1,j=1}^{M,N} \left( \left( \mathrm{D} \mathbf{x} \right)_{i,j,1}^p + \left( \mathrm{D} \mathbf{x} \right)_{i,j,2}^p \right)^{1/p}" />

where <img src="https://latex.codecogs.com/svg.latex?\inline&space;\large&space;\mathrm{D}:&space;\mathbb{R}^{M&space;\times&space;N}&space;\rightarrow&space;\mathbb{R}^{M&space;\times&space;N&space;\times&space;2}" title="\large \mathrm{D}: \mathbb{R}^{M \times N} \rightarrow \mathbb{R}^{M \times N \times 2}" /> is the discrete gradient operator. That is, <img src="https://latex.codecogs.com/svg.latex?\inline&space;\large&space;\left\Vert&space;\cdot&space;\right\Vert_{\text{TV}}" title="\large \left\Vert \cdot \right\Vert_{\text{TV}}" /> is the <img src="https://latex.codecogs.com/svg.latex?\inline&space;\large&space;\ell_1" title="\large \ell_1" />-norm  of the <img src="https://latex.codecogs.com/svg.latex?\inline&space;\large&space;p" title="\large p" />-norm of the pixelwise image gradients [1, pp. 168]. When <img src="https://latex.codecogs.com/svg.latex?\inline&space;\large&space;p=1" title="\large p=1" /> it is called the anisotropic TV, whereas when <img src="https://latex.codecogs.com/svg.latex?\inline&space;\large&space;p=2" title="\large p=2" /> it is called the isotropic TV, which is used here.

The code solves the problem using the over-relaxed Chambolle-Pock algorithm [2]

## References
- A. Chambolle and T. Pock, “An introduction to continuous optimization for imaging,” _Acta Numer._, vol. 25, pp. 161–319, 2016.

<!--stackedit_data:
eyJoaXN0b3J5IjpbLTE4NzE0MzE1ODQsLTgwOTUwMDM1MCwtMT
YwODcyOTY2NSwyMzQ0NzIxNDcsMjA0OTE5NTkxMiwxNTEyODYx
NjU1LC0yMTQ3MzU1OCwtMTkwODYxNDcxMiwtMjY2Mjc0OTkxLD
EzNjk5OTg1NzUsLTE3Njk2MTEzNzksLTE3NTc4NTkwOTAsLTI0
NjYxNzc4Ml19
-->