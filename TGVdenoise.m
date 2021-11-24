%%
% Denoising/smoothing a given image y with the second order 
% total generalized variation (TGV)
%
% The iterative algorithm converges to the unique image u 
% (and the auxiliary vectorial image v) minimizing
%
% ||u-y||_2^2/2 + lambda0*||Jv||_2,1 + lambda1*||Du-v||_2,1 
%
% where D maps an image to its gradient field (i.e. R^{HxW} -> R^{HxWx2})
% and J maps a vector field to its Jacobian (i.e. R^{HxWx2} -> R^{HxWx4})
%		
% This code is adapted from L. Condat's TGVdenoise.m
% The original code is available from
% https://lcondat.github.io/software.html
% 
% Changes include:
% 1. The problem is reformulated
% 2. opJ is changed to obey Dirichlet boundary condition
% 3. opJadj is changed to obey Dirichlet boundary condition accordingly
%
% see
% <An introduction to continuous optimization for imaging> pp. 245
% for the problem formulation
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

	Nbiter= 600;	% number of iterations
	lambda0 = 1/4; 	% regularization parameter
	lambda1 = 1/9;	% regularization parameter
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
    
    [xsol, primal_cost, dual_cost] = TGVdenoising_TD(y,lambda0,lambda1,tau,Nbiter);
	subplot(223);
	imshow(xsol);
    title('TGV denoised image');
    imwrite(xsol,'images\TGVdenoised_gray.png');
    
    subplot(224);
    plot(primal_cost);
    xlabel('iteration');
    grid on;
    hold on;
    plot(dual_cost)
    title('Primal and dual cost');

    fprintf('noisy image: RSNR = %.4f dB, SSIM = %.4f\n',calcRSNR(y,x0),ssim(y,x0));
    fprintf('TGV denoised image: RSNR = %.4f dB, SSIM = %.4f\n',calcRSNR(xsol,x0),ssim(xsol,x0));
end


function [u, primal_cost, dual_cost] = TGVdenoising_TD(y,lambda0,lambda1,tau,Nbiter)
    
    rho = 1.99;		% relaxation parameter, in [1,2)
	sigma = 1/tau/72; % the dual step size
	[H,W] = size(y);
    
    opD = @(u) cat(3,[diff(u,1,1);zeros(1,W)],[diff(u,1,2) zeros(H,1)]);    % grad with Neumann boundary condition
	opDadj = @(p) [-p(1,:,1);-diff(p(1:end-1,:,1),1,1);p(end-1,:,1)]+...    % -div with Dirichlet boundary condition
		[-p(:,1,2) -diff(p(:,1:end-1,2),1,2) p(:,end-1,2)];	
    opJ = @(v) cat(3,...                                                    % grad with Neumann boundary condition applied to a vectorial image: R^(MxNx2) -> R^(MxNx4)
		[diff(v(:,:,1),1,1);zeros(1,W)],[diff(v(:,:,1),1,2) zeros(H,1)],...
		[diff(v(:,:,2),1,1);zeros(1,W)],[diff(v(:,:,2),1,2) zeros(H,1)]);
	opJadj = @(q) cat(3,...
		[-q(1,:,1);-diff(q(1:end-1,:,1),1,1);q(end-1,:,1)]+[-q(:,1,2) -diff(q(:,1:end-1,2),1,2) q(:,end-1,2)],...  % -div with Dirichlet boundary condition applied to a vectorial image: R^(MxNx4) -> R^(MxNx2)
		[-q(1,:,3);-diff(q(1:end-1,:,3),1,1);q(end-1,:,3)]+[-q(:,1,4) -diff(q(:,1:end-1,4),1,2) q(:,end-1,4)]);
    prox_tau_g = @(u) (u+tau*y)/(1+tau);                                    % prox of tau*||.-y||_2^2/2
    prox_gamma_f = @(u,gamma) max(1-gamma./sqrt(sum(u.^2,3)),0).*u;         % prox of gamma*||.||_2
    
    u2 = y;                 % Initialization of the solution
	v2 = zeros([H,W,2]);	% Initialization of the vector field v
    p2 = zeros([H,W,2]);	% Initialization of the dual solution p
	q2 = zeros([H,W,4]);	% Initialization of the dual solution q
    
    primal_cost = NaN(1, Nbiter);
    dual_cost = NaN(1, Nbiter);

    for iter = 1:Nbiter
        
        % primal updates
        u = prox_tau_g(u2 - tau*opDadj(p2));
        v = v2 - tau*(-p2 + opJadj(q2));
        
        % dual updates using Moreau proximal decomposition theorem
        p = p2 + sigma*(opD(2*u-u2)-(2*v-v2));
        p = p - sigma*prox_gamma_f(p/sigma,lambda1/sigma);
        q = q2 + sigma*opJ(2*v-v2);
        q = q - sigma*prox_gamma_f(q/sigma,lambda0/sigma);
        
        % over-relaxations
        u2 = u2 + rho*(u-u2);
		v2 = v2 + rho*(v-v2);
		p2 = p2 + rho*(p-p2);
        q2 = q2 + rho*(q-q2);
        
        % evaluate cost values
        primal_cost(iter) = sum((u-y).^2,'all')/2+...
				lambda0*sum(sum(sqrt(sum(opJ(v).^2,3))))+...
				lambda1*sum(sum(sqrt(sum((opD(u)-v).^2,3))));
        dual_cost(iter) = 0.5*sum((cat(3,y,zeros([H,W,2]))).^2,'all') ...
                        - 0.5*sum((cat(3,y,zeros([H,W,2]))-cat(3,opDadj(p),opJadj(q)-p)).^2,'all');
        if mod(iter,25)==0
			fprintf('nb iter:%4d  %f  %f\n',iter,primal_cost(iter),dual_cost(iter));
        end
        
    end
end
