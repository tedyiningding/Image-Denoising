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
    
    x0 = double(imread('images\parrotgray.png'))/255;   % Initial image
	figure(1);
	imshow(x0);
    title('clean image')
    
	rng(0);
	y = x0+randn(size(x0))*0.1; % white Gaussian noise added to the image
	figure(2);
	imshow(y);
    title(['noisy image with RSNR = ',num2str(calcRSNR(y,x0))]);
	imwrite(y,'images\noisy.png');
    
    xsol = TGVdenoising_TD(y,lambda0,lambda1,tau,Nbiter);
	figure(3);
	imshow(xsol);
    title(['denoised image with RSNR = ',num2str(calcRSNR(xsol,x0))]);
    imwrite(xsol,'images\TGVdenoised.png');
    
end


function u = TGVdenoising_TD(y,lambda0,lambda1,tau,Nbiter)
    
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
    
	cy = sum(sum(y.^2))/2;
	primalcostlowerbound = 0;
    
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
        
        if mod(iter,25)==0
			primalcost = norm(u-y,'fro')^2/2+...
				lambda0*sum(sum(sqrt(sum(opJ(v).^2,3))))+...
				lambda1*sum(sum(sqrt(sum((opD(u)-v).^2,3))));		
			dualcost = cy-sum(sum((y-opDadj(opJadj(q))).^2))/2;
			tmp = max(max(sqrt(sum(opJadj(q).^2,3))));
			  % the value tmp is to check feasibility: the value will be
			  % <= lambda1 only at convergence. Since u is not feasible,
			  % the dual cost is not reliable: the gap=primalcost-dualcost
			  % can be <0 and cannot be used as stopping criterion.
			q3=q/max(tmp/lambda0,1);
			  % q3 is a scaled version of u, which is feasible.
			  % so, its dual cost is a valid, but very rough lower bound
			  % of the primal cost. 
			dualcost2 = cy-sum(sum((y-opDadj(opJadj(q3))).^2))/2;
			  % we display the best value of dualcost2 computed so far.
			primalcostlowerbound = max(primalcostlowerbound,dualcost2);
			  % The gap between primalcost and primalcostlowerbound tends 
			  % to zero but does not tell much about convergence, since x 
			  % converges much faster than u3.
			fprintf('nb iter:%4d  %f  %f  %f  %f\n',iter,primalcost,...
				dualcost,primalcostlowerbound,tmp);
			figure(3);
			imshow(u);
        end
        
    end
end
