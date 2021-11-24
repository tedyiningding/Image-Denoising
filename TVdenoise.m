%%
% Denoising/smoothing a given image y with the isotropic total variation.
%
% The iterative algorithm converges to the unique image x minimizing 
%
% ||x-y||_2^2/2 + lambda*TV(x)
%     g(x)      +     f(x)
% 
% TV(x) = ||Dx||_1,2, where D maps an image to its gradient field.
%
% The code is adapted from L. Condat's TVdenoise.m
% The original code is available from
% https://lcondat.github.io/software.html
% 
% Changes include:
% 1. opDadj (i.e. -div) is changed to obey Dirichlet boundary condition
% 2. Use prox_sigma_g together with Moreau proximal decomposition theorem, 
%    instead of directly applying prox_sigma_g_conj
%
% see
% <An introduction to continuous optimization for imaging> pp. 52, 58, 59
% for the algorithm
% and
% <A primal-dual splitting method for convex optimization 
% involving Lipschitzian, proximable and linear composite terms>
% pp. 4 Algorithm 3.1 for the over-relaxations
%
% Ted (Yining) Ding, PhD student in Robotics and Autonomous Systems
% Edinburgh Centre for Robotics
% yd2007@hw.ac.uk
% 16/07/2021
%%
function main
    
    addpath('utils\')

	Nbiter= 400;	% number of iterations
	lambda = 0.1; 	% regularization parameter
	tau = 0.01;		% proximal parameter > 0; influences the convergence speed (i.e. primal step size)
	
    figure;

	x0 = double(imread('images\gray.png'))/255;   % Initial image
	subplot(221);
	imshow(x0);
    title('clean image')
    
	rng(0);
	y = x0+randn(size(x0))*0.1; % white Gaussian noise added to the image
	subplot(222);
	imshow(y);
    title('noisy image');
    imwrite(y,'images\noisy_gray.png');
    
	[xsol, primal_cost, dual_cost] = TVdenoising(y,lambda,tau,Nbiter);
	subplot(223);
	imshow(xsol);
    title('TV denoised image');
    imwrite(xsol,'images\TVdenoised_gray.png');
    
    subplot(224);
    plot(primal_cost);
    xlabel('iteration');
    grid on;
    hold on;
    plot(dual_cost)
    title('Primal and dual cost');

    fprintf('noisy image: RSNR = %.4f dB, SSIM = %.4f\n',calcRSNR(y,x0),ssim(y,x0));
    fprintf('TV denoised image: RSNR = %.4f dB, SSIM = %.4f\n',calcRSNR(xsol,x0),ssim(xsol,x0));
end


function [x, primal_cost, dual_cost] = TVdenoising(y,lambda,tau,Nbiter)
	
	rho = 1.99;         % relaxation parameter, in [1,2)
	sigma = 1/tau/8;    % the dual step size (N.B. ||D||^2 = 8)

	opD = @(x) cat(3,[diff(x,1,1);zeros(1,size(x,2))],[diff(x,1,2) zeros(size(x,1),1)]);            % grad with Neumann boundary condition
    opDadj = @(u) -[u(1,:,1);u(2:end-1,:,1)-u(1:end-2,:,1);-u(end-1,:,1)]-[u(:,1,2) u(:,2:end-1,2)-u(:,1:end-2,2) -u(:,end-1,2)];     % -div with Dirichlet boundary condition
	prox_tau_g = @(x) (x+tau*y)/(1+tau);                                            % prox of tau*||.-y||_2^2/2
    prox_gamma_f = @(u,gamma) max(1-gamma./sqrt(sum(u.^2,3)),0).*u;                 % prox of gamma*||.||_2
	
	x2 = y; 		% Initialization of the solution
    u2 = zeros([size(y),2]);
    
	cy = sum(sum(y.^2))/2;
	
    primal_cost = NaN(1, Nbiter);
    dual_cost = NaN(1, Nbiter);
		
    for iter = 1:Nbiter
        
        % primal update
		x = prox_tau_g(x2 - tau*opDadj(u2));
        
        % dual update using Moreau proximal decomposition theorem
        u = u2 + sigma*opD(2*x - x2);
        u = u - sigma*prox_gamma_f(u/sigma,lambda/sigma);
        
		% over-relaxations
        x2 = x2 + rho*(x-x2);
		u2 = u2 + rho*(u-u2);
        
        % evaluate cost values
        primal_cost(iter) = sum((x-y).^2,'all')/2 + lambda*sum(sum(sqrt(sum(opD(x).^2,3))));    % ||x-y||_2^2/2 + lambda.TV(x)
        dual_cost(iter) = cy - sum(sum((y-opDadj(u)).^2))/2;                                    % ||y||_2^2/2 - ||y+div(u)||^2/2, which is consistent with <An introduction to continuous optimization for imaging> pp 18
        
        % to monitor convergence
        if mod(iter,25)==0
            fprintf('nb iter:%4d  %f  %f\n',iter,primal_cost(iter),dual_cost(iter));
        end
        
    end
end
