
# Image-Denoising
A collection of MATLAB scripts and functions which achieves image denoising based on variational methods

Some codes are adapted from [Laurent Condat's work](https://lcondat.github.io/software.html).

## TV Denoising (TVdenoise.m)
TV denoising finds the unique solution <img src="https://latex.codecogs.com/svg.latex?\inline&space;\large&space;\mathbf{x}^\star&space;\in&space;\mathbb{R}^{M&space;\times&space;N}" title="\large \mathbf{x}^\star \in \mathbb{R}^{M \times N}" /> which minimises

<img src="https://latex.codecogs.com/svg.latex?\large&space;\frac{1}{2}&space;\left\Vert&space;\mathbf{x}-\mathbf{y}&space;\right\Vert&space;^2&space;&plus;&space;\lambda&space;\left\Vert&space;\mathbf{x}&space;\right\Vert_{\text{TV}}" title="\large \frac{1}{2} \left\Vert \mathbf{x}-\mathbf{y} \right\Vert ^2 + \lambda \left\Vert \mathbf{x} \right\Vert_{\text{TV}}" />

where <img src="https://latex.codecogs.com/svg.latex?\inline&space;\large&space;\mathbf{y}&space;\in&space;\mathbb{R}^{M&space;\times&space;N}" title="\large \mathbf{y} \in \mathbb{R}^{M \times N}" /> is the observed image, <img src="https://latex.codecogs.com/svg.latex?\inline&space;\large&space;\lambda" title="\large \lambda" /> is a regularisation parameter that balances the two terms, and <img src="https://latex.codecogs.com/svg.latex?\inline&space;\large&space;\left\Vert&space;\mathbf{x}&space;\right\Vert_{\text{TV}}" title="\large \left\Vert \mathbf{x} \right\Vert_{\text{TV}}" /> is the total variation which is defined below

<img src="https://latex.codecogs.com/svg.latex?\large&space;\left\Vert&space;\mathbf{x}&space;\right\Vert_{\text{TV}}&space;=&space;\left\Vert&space;\mathrm{D}&space;\mathbf{x}&space;\right\Vert_{p,1}&space;=&space;\sum_{i=1,j=1}^{M,N}&space;\left|&space;\left(&space;\mathrm{D}&space;\mathbf{x}&space;\right)_{i,j}&space;\right|_p&space;=&space;\sum_{i=1,j=1}^{M,N}&space;\left(&space;\left(&space;\mathrm{D}&space;\mathbf{x}&space;\right)_{i,j,1}^p&space;&plus;&space;\left(&space;\mathrm{D}&space;\mathbf{x}&space;\right)_{i,j,2}^p&space;\right)^{1/p}" title="\large \left\Vert \mathbf{x} \right\Vert_{\text{TV}} = \left\Vert \mathrm{D} \mathbf{x} \right\Vert_{p,1} = \sum_{i=1,j=1}^{M,N} \left| \left( \mathrm{D} \mathbf{x} \right)_{i,j} \right|_p = \sum_{i=1,j=1}^{M,N} \left( \left( \mathrm{D} \mathbf{x} \right)_{i,j,1}^p + \left( \mathrm{D} \mathbf{x} \right)_{i,j,2}^p \right)^{1/p}" />


<!--stackedit_data:
eyJoaXN0b3J5IjpbMzM1NTQxODkxLDE1MTI4NjE2NTUsLTIxND
czNTU4LC0xOTA4NjE0NzEyLC0yNjYyNzQ5OTEsMTM2OTk5ODU3
NSwtMTc2OTYxMTM3OSwtMTc1Nzg1OTA5MCwtMjQ2NjE3NzgyXX
0=
-->