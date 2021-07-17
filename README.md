
# Image-Denoising
A collection of MATLAB scripts and functions which achieves image denoising based on variational methods

Some codes are adapted from [Laurent Condat's work](https://lcondat.github.io/software.html).

## TV Denoising
TV denoising finds the unique solution <img src="https://latex.codecogs.com/svg.latex?\inline&space;\large&space;\mathbf{x}^\star&space;\in&space;\mathbb{R}^{M&space;\times&space;N}" title="\large \mathbf{x}^\star \in \mathbb{R}^{M \times N}" /> which minimises

<img src="https://latex.codecogs.com/svg.latex?\large&space;\frac{1}{2}&space;\left\Vert&space;\mathbf{x}-\mathbf{y}&space;\right\Vert&space;^2&space;&plus;&space;\lambda&space;\left\Vert&space;\mathbf{x}&space;\right\Vert_{\text{TV}}" title="\large \frac{1}{2} \left\Vert \mathbf{x}-\mathbf{y} \right\Vert ^2 + \lambda \left\Vert \mathbf{x} \right\Vert_{\text{TV}}" />

where <img src="https://latex.codecogs.com/svg.latex?\inline&space;\large&space;\mathbf{y}&space;\in&space;\mathbb{R}^{M&space;\times&space;N}" title="\large \mathbf{y} \in \mathbb{R}^{M \times N}" /> is the observed image, <img src="https://latex.codecogs.com/svg.latex?\inline&space;\large&space;\lambda" title="\large \lambda" /> is a regularisation parameter that balances the two terms, and <img src="https://latex.codecogs.com/svg.latex?\inline&space;\large&space;\left\Vert&space;\mathbf{x}&space;\right\Vert_{\text{TV}}" title="\large \left\Vert \mathbf{x} \right\Vert_{\text{TV}}" /> is the isotropic total variation which is defined below



$$\frac{1}{2} \left\Vert \mathbf{x}-\mathbf{y} \right\Vert ^2 + \lambda \left\Vert \mathbf{x} \right\Vert_{\text{TV}}$$
where

<!--stackedit_data:
eyJoaXN0b3J5IjpbLTE5MDg2MTQ3MTIsLTI2NjI3NDk5MSwxMz
Y5OTk4NTc1LC0xNzY5NjExMzc5LC0xNzU3ODU5MDkwLC0yNDY2
MTc3ODJdfQ==
-->