# Image-Denoising
A collection of MATLAB scripts and functions which achieves image denoising based on variational methods

Some codes are adapted from [Laurent Condat's work](https://lcondat.github.io/software.html).

## TV Denoising
TV denoising finds the unique solution $\mathbf{x}^\star \in \mathbb{R}^{M \times N}$ which minimises
$$\frac{1}{2} \left\Vert \mathbf{x}-\mathbf{y} \right\Vert ^2 + \lambda \left\Vert \mathbf{x} \right\Vert_{\text{TV}}$$
where

<!--stackedit_data:
eyJoaXN0b3J5IjpbMTM2OTk5ODU3NSwtMTc2OTYxMTM3OSwtMT
c1Nzg1OTA5MCwtMjQ2NjE3NzgyXX0=
-->