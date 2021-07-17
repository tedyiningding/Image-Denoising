% This code is adapted from L. Condat's TGVdenoise.m by Ted Ding.
% 30/06/2021
% The original code is available from
% https://lcondat.github.io/software.html
% 
% Changes include:
% 1. opJ is changed to obey Dirichlet boundary condition
% 2. opJadj is changed to obey Dirichlet boundary condition accordingly
%
% Denoising/smoothing a given image y with the second order 
% total generalized variation (TGV), defined in 
% K. Bredies, K. Kunisch, and T. Pock, "Total generalized variation,"
% SIAM J. Imaging Sci., 3(3), 492-526, 2010.
%
% The iterative algorithm converges to the unique image u 
% (and the auxiliary vectorial image v) minimizing the following defined in
% <An introduction to continuous optimization for imaging> pp. 245
%
% lambda1.||Du-v||_2,1 + lambda0.||Jv||_2,1 + ||u-y||_2^2/2  
% 
%
% where D maps an image to its gradient field (i.e. R^{MxN} -> R^{MxNx2})
% and J maps a vector field to its Jacobian (i.e. R^{MxNx2} -> R^{MxNx4})
%		
% The over-relaxed Chambolle-Pock algorithm is described in
% L. Condat, "A primal-dual splitting method for convex optimization 
% involving Lipschitzian, proximable and linear composite terms", 
% J. Optimization Theory and Applications, vol. 158, no. 2, 
% pp. 460-479, 2013.
%
% Code written by Laurent Condat, CNRS research fellow in the
% Dept. of Images and Signals of GIPSA-lab, Univ. Grenoble Alpes, 
% Grenoble, France.
%
% Version 1.0, Oct. 12, 2016


function main

	Nbiter= 600;	% number of iterations
	lambda0 = 0.1; 	% regularization parameter
	lambda1 = 0.15;	% regularization parameter
	tau = 0.01;		% proximal parameter >0; influences the
		% convergence speed

	y  = double(imread('parrotgray.png'))/255;   % Initial image
	figure(1);
	imshow(y);
	rng(0);
	y = y+randn(size(y))*0.1; % white Gaussian noise added to the image
	figure(2);
	imshow(y);
	x = TGVdenoising_TD(y,lambda0,lambda1,tau,Nbiter);
	figure(3);
	imshow(x);
	imwrite(y,'noisy.png');
	imwrite(x,'TGVdenoised_TD.png');
end


function u = TGVdenoising_TD(y,lambda0,lambda1,tau,Nbiter)
    
    rho = 1.99;		% relaxation parameter, in [1,2)
	sigma = 1/tau/72; % proximal parameter
	[H,W]=size(y);
    
    % step sizes - primal (tau) and dual (sigma)
    tau_u = 1;
    tau_v = 1;
    sigma_p = 1;
    sigma_q = 1;
    
    opD = @(u) cat(3,[diff(u,1,1);zeros(1,W)],[diff(u,1,2) zeros(H,1)]);    % TD: grad with Neumann boundary condition
	opDadj = @(p) [-p(1,:,1);-diff(p(1:end-1,:,1),1,1);p(end-1,:,1)]+...    % TD: -div with Dirichlet boundary condition
		[-p(:,1,2) -diff(p(:,1:end-1,2),1,2) p(:,end-1,2)];	
    opJ = @(v) cat(3,...                                                    % TD: grad with Neumann boundary condition applied to a vectorial image: R^(MxNx2) -> R^(MxNx4)
		[diff(v(:,:,1),1,1);zeros(1,W)],[diff(v(:,:,1),1,2) zeros(H,1)],...
		[diff(v(:,:,2),1,1);zeros(1,W)],[diff(v(:,:,2),1,2) zeros(H,1)]);
	opJadj = @(q) cat(3,...
		[-q(1,:,1);-diff(q(1:end-1,:,1),1,1);q(end-1,:,1)]+[-q(:,1,2) -diff(q(:,1:end-1,2),1,2) q(:,end-1,2)],...  % TD: -div with Dirichlet boundary condition applied to a vectorial image: R^(MxNx4) -> R^(MxNx2)
		[-q(1,:,3);-diff(q(1:end-1,:,3),1,1);q(end-1,:,3)]+[-q(:,1,4) -diff(q(:,1:end-1,4),1,2) q(:,end-1,4)]);
    prox_tau_u_g = @(u) (u+tau*y)/(1+tau); 
%     prox_gamma_f = @(u,gamma) max(1-gamma./sqrt(sum(u.^2,3)),0).*u;                     % TD: prox of gamma*||.||_2
    prox_sigma_f_conj = @(u) bsxfun(@rdivide,u,max(sqrt(sum(u.^2,3))/...    % TD: projection onto l2-ball centred at 0 with radius lambda (i.e. prox of indicator function of l2-ball centred at 0 with radius lambda2)
		1,1));
    
    u2 = y; 		% Initialization of the solution
	v2 = zeros([H,W,2]);	% Initialization of the vector field v
    p2 = zeros([H,W,2]);	% Initialization of the dual solution p
	q2 = zeros([H,W,4]);	% Initialization of the dual solution q
	cy = sum(sum(y.^2))/2;
	primalcostlowerbound = 0;
    
    for iter = 1:Nbiter
        % primal updates
        u = prox_tau_u_g(u2-tau*lambda1*opDadj(p2));
        v = v2-tau*(-lambda1*p2+lambda0*opJadj(q2));
        % dual updates
%         p = p2+sigma*lambda1*(opD(2*u-u2)-v);
%         p = p - sigma*prox_gamma_f(p/sigma,1/sigma);
%         q = q2+sigma*lambda0*opJ(2*v-v2);
%         q = q - sigma*prox_gamma_f(q/sigma,1/sigma);
        p = prox_sigma_f_conj(p2+sigma*lambda1*(opD(2*u-u2)-(2*v-v2)));
        q = prox_sigma_f_conj(q2+sigma*lambda0*opJ(2*v-v2));
        % over-relaxation
        u2 = u2+rho*(u-u2);
		v2 = v2+rho*(v-v2);
		p2 = p2+rho*(p-p2);
        q2 = q2+rho*(q-q2);
        
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
