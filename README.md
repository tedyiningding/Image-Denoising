
# Image-Denoising
A collection of MATLAB scripts and functions which achieves image denoising based on variational methods

Some codes are adapted from [Laurent Condat's work](https://lcondat.github.io/software.html).

## TV Denoising
Total Variation (TV) denoising (more specifically, we mean the ROF model [1]) finds the unique solution <img src="https://latex.codecogs.com/svg.latex?\inline&space;\large&space;\mathbf{x}^\star&space;\in&space;\mathbb{R}^{M&space;\times&space;N}" title="\large \mathbf{x}^\star \in \mathbb{R}^{M \times N}" /> which minimises

<img src="https://latex.codecogs.com/svg.latex?\large&space;\frac{1}{2}&space;\left\Vert&space;\mathbf{x}-\mathbf{y}&space;\right\Vert&space;^2&space;&plus;&space;\lambda&space;\left\Vert&space;\mathbf{x}&space;\right\Vert_{\text{TV}}" title="\large \frac{1}{2} \left\Vert \mathbf{x}-\mathbf{y} \right\Vert ^2 + \lambda \left\Vert \mathbf{x} \right\Vert_{\text{TV}}" />

where <img src="https://latex.codecogs.com/svg.latex?\inline&space;\large&space;\mathbf{y}&space;\in&space;\mathbb{R}^{M&space;\times&space;N}" title="\large \mathbf{y} \in \mathbb{R}^{M \times N}" /> is the observed image, <img src="https://latex.codecogs.com/svg.latex?\inline&space;\large&space;\lambda" title="\large \lambda" /> is a regularisation parameter that balances the two terms, and <img src="https://latex.codecogs.com/svg.latex?\inline&space;\large&space;\left\Vert&space;\cdot&space;\right\Vert_{\text{TV}}" title="\large \left\Vert \cdot \right\Vert_{\text{TV}}" /> is the total variation of an image which is defined by

<img src="https://latex.codecogs.com/svg.latex?\large&space;\left\Vert&space;\mathbf{x}&space;\right\Vert_{\text{TV}}&space;=&space;\left\Vert&space;\mathrm{D}&space;\mathbf{x}&space;\right\Vert_{p,1}&space;=&space;\sum_{i=1,j=1}^{M,N}&space;\left|&space;\left(&space;\mathrm{D}&space;\mathbf{x}&space;\right)_{i,j}&space;\right|_p&space;=&space;\sum_{i=1,j=1}^{M,N}&space;\left(&space;\left(&space;\mathrm{D}&space;\mathbf{x}&space;\right)_{i,j,1}^p&space;&plus;&space;\left(&space;\mathrm{D}&space;\mathbf{x}&space;\right)_{i,j,2}^p&space;\right)^{1/p}" title="\large \left\Vert \mathbf{x} \right\Vert_{\text{TV}} = \left\Vert \mathrm{D} \mathbf{x} \right\Vert_{p,1} = \sum_{i=1,j=1}^{M,N} \left| \left( \mathrm{D} \mathbf{x} \right)_{i,j} \right|_p = \sum_{i=1,j=1}^{M,N} \left( \left( \mathrm{D} \mathbf{x} \right)_{i,j,1}^p + \left( \mathrm{D} \mathbf{x} \right)_{i,j,2}^p \right)^{1/p}" />

where <img src="https://latex.codecogs.com/svg.latex?\inline&space;\large&space;\mathrm{D}:&space;\mathbb{R}^{M&space;\times&space;N}&space;\rightarrow&space;\mathbb{R}^{M&space;\times&space;N&space;\times&space;2}" title="\large \mathrm{D}: \mathbb{R}^{M \times N} \rightarrow \mathbb{R}^{M \times N \times 2}" /> is the discrete gradient operator. That is, <img src="https://latex.codecogs.com/svg.latex?\inline&space;\large&space;\left\Vert&space;\cdot&space;\right\Vert_{\text{TV}}" title="\large \left\Vert \cdot \right\Vert_{\text{TV}}" /> is the <img src="https://latex.codecogs.com/svg.latex?\inline&space;\large&space;\ell_1" title="\large \ell_1" />-norm  of the <img src="https://latex.codecogs.com/svg.latex?\inline&space;\large&space;p" title="\large p" />-norm of the pixelwise image gradients [2, pp. 168]. When <img src="https://latex.codecogs.com/svg.latex?\inline&space;\large&space;p=1" title="\large p=1" /> it is called the anisotropic TV, whereas when <img src="https://latex.codecogs.com/svg.latex?\inline&space;\large&space;p=2" title="\large p=2" /> it is called the isotropic TV. The code uses the latter.

The [code](https://github.com/tedyiningding/Image-Denoising/blob/main/TVdenoise.m) solves the problem using the over-relaxed Chambolle-Pock algorithm [3, Algorithm 3.1] after obtaining a saddle-point problem [2, Example 5.6]

## TGV denoising
TV regularisation only promotes piecewise constant therefore the result could suffer from staircasing artefacts (as will be seen from the denoised images shown below). To combat this, Total Generalised Variation (TGV) was proposed in [4]. The second order TGV not only promotes piecewise constant, but also piecewise affine. It finds the unique solution <img src="https://latex.codecogs.com/svg.latex?\inline&space;\large&space;\mathbf{u}^\star&space;\in&space;\mathbb{R}^{M&space;\times&space;N}" title="\large \mathbf{u}^\star \in \mathbb{R}^{M \times N}" /> (and an auxiliary variable <img src="https://latex.codecogs.com/svg.latex?\inline&space;\large&space;\mathbf{v}^\star&space;\in&space;\mathbb{R}^{M&space;\times&space;N&space;\times&space;2}" title="\large \mathbf{v}^\star \in \mathbb{R}^{M \times N \times 2}" />) which minimises



## Results
The clear image | The noisy image
:-:|:-:
![clear](https://github.com/tedyiningding/Image-Denoising/blob/main/images/gray.png?raw=true) | ![noisy](https://github.com/tedyiningding/Image-Denoising/blob/main/images/noisy_gray.png?raw=true)
TV denoised image | TGV denoised image
![clear](https://github.com/tedyiningding/Image-Denoising/blob/main/images/TVdenoised_gray.png?raw=true) | ![noisy](https://github.com/tedyiningding/Image-Denoising/blob/main/images/TGVdenoised_gray.png?raw=true)

The table below shows some quantitative results including the reconstructed signal-to-noise ratio (RSNR) and the structural similarity index measure (SSIM).
|  Images | RSNR (dB) | SSIM |
| - | - | - |
| Noisy        | 14.1727| 0.1949
| TV denoised  | 25.3879| 0.8702
| TGV denoised | 25.3985| 0.8859

The TGV method is only marginally higher in RSNR and SSIM but there is no obvious staircasing artefact (see the background).

## References
- [1] L. I. Rudin, S. Osher, and E. Fatemi, “Nonlinear total variation based noise removal algorithms,” _Physica D_, vol. 60, no. 1–4, pp. 259–268, 1992.
- [2] A. Chambolle and T. Pock, “An introduction to continuous optimization for imaging,” _Acta Numer._, vol. 25, pp. 161–319, 2016.
- [3] L. Condat, “A primal–dual splitting method for convex optimization involving lipschitzian, proximable and linear composite terms,” _J. Optim. Theory Appl._, vol. 158, no. 2, pp. 460–479, 2013.
- [4] K. Bredies, K. Kunisch, and T. Pock, “Total Generalized Variation,” _SIAM J. Imaging Sci._, vol. 3, no. 3, pp. 492–526, 2010.

<!--stackedit_data:
eyJoaXN0b3J5IjpbLTE1OTUyNTM5NDIsMjkyODkwMTg2LC0zMz
cxMjAxNTEsLTYxMzY2MDg4LC00MDMxMTUyOTcsLTIwNDc0OTQz
OTEsNjA4NDM0OTEwLC0xODcxNDMxNTg0LC04MDk1MDAzNTAsLT
E2MDg3Mjk2NjUsMjM0NDcyMTQ3LDIwNDkxOTU5MTIsMTUxMjg2
MTY1NSwtMjE0NzM1NTgsLTE5MDg2MTQ3MTIsLTI2NjI3NDk5MS
wxMzY5OTk4NTc1LC0xNzY5NjExMzc5LC0xNzU3ODU5MDkwLC0y
NDY2MTc3ODJdfQ==
-->