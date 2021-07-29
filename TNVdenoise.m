% Denoising/smoothing a given color image y with the 
% isotropic total nuclear variation.
%
% The iterative algorithm converges to the unique image x minimizing 
%
% ||x-y||_2^2/2 + lambda.TNV(x)
%
% TV(x)=||Dx||_1,*, where D maps a color image to its Jacobian field
% and * is the nuclear norm (of the 3x2 Jacobian matrix at each pixel)
%
% This penaly was called the total nuclear variation in
% * K.M. Holt, Total nuclear variation and Jacobian extensions
% of total variation for vector fields, IEEE Trans. Image Proc.,
% vol. 23, pp. 3975–3989, 2014
% It has also been studied with other names in
% * S. Lefkimmiatis, A. Roussos, M. Unser, and P. Maragos, Convex 
% generalizations of total variation based on the structure tensor
% with applications to inverse problems, in Scale Space and 
% Variational Methods in Computer Vision, Lecture Notes in Comput.
% Sci. 7893, Springer, Berlin, 2013, pp. 48–60.
% * G Chierchia, N Pustelnik, B Pesquet-Popescu, JC Pesquet, 
% "A nonlocal structure tensor-based approach for multicomponent
% image recovery problems", IEEE Trans. Image Proc., 23 (12), 
% pp. 5531-5544, 2014.
% * J. Duran, M. Moeller, C. Sbert, and D. Cremers, "Collaborative
% Total Variation: A General Framework for Vectorial TV Models", 
% SIAM J. Imaging Sciences, Vol. 9, No. 1, pp. 116–151, 2016.
%
%
% The over-relaxed Chambolle-Pock algorithm used here is described
% in L. Condat, "A primal-dual splitting method for convex
% optimization involving Lipschitzian, proximable and linear
% composite terms", J. Optimization Theory and Applications, 
% vol. 158, no. 2, pp. 460-479, 2013.
%
% Code written by Laurent Condat, CNRS research fellow in the
% Dept. of Images and Signals of GIPSA-lab, Univ. Grenoble Alpes, 
% Grenoble, France.
%
% Version 1.0, Jul. 12, 2018


function main

	Nbiter= 1000;	% number of iterations
	lambda = 0.12; 	% regularization parameter
	tau = 0.005;		% proximal parameter >0; influences the
		% convergence speed
				
	y  = double(imread('images\colour.tif'))/255;   % Initial image
	figure(1);
	imshow(y);
	rng(0);
	y = y+randn(size(y))*0.1; % white Gaussian noise added to the image
	figure(2);
	imshow(y);
	x = TNVdenoising(y,lambda,tau,Nbiter);
	figure(3);
	imshow(x);
	imwrite(y,'images\noisy_colour.png');
	imwrite(x,'images\TNVdenoised_colour.png');
end


function x = TNVdenoising(y,lambda,tau,Nbiter)
	
	rho = 1.99;		% relaxation parameter, in [1,2)
	sigma = 1/tau/8; % proximal parameter
	[H,W,C]=size(y);

	opD = @(x) cat(4,[diff(x,1,1);zeros(1,W,C)],[diff(x,1,2) zeros(H,1,C)]);                % R^(MxNx3) -> R^(MxNx3x2)
	opDadj = @(u) -[u(1,:,:,1);u(2:end-1,:,:,1)-u(1:end-2,:,:,1);-u(end-1,:,:,1)]-[u(:,1,:,2) u(:,2:end-1,:,2)-u(:,1:end-2,:,2) -u(:,end-1,:,2)];     % R^(MxNx3x2) -> R^(MxNx3)
    prox_tau_f = @(x) (x+tau*y)/(1+tau);
    
	x2 = y; 		% Initialization of the solution
	u2 = zeros([size(y) 2]); % Initialization of the dual solution
	cy = sum(sum(sum(y.^2)))/2;
	primalcostlowerbound = 0;
		
	for iter = 1:Nbiter
		x = prox_tau_f(x2-tau*opDadj(u2));
        u = prox_sigma_g_conj(u2+sigma*opD(2*x-x2),lambda);
		x2 = x2+rho*(x-x2);
		u2 = u2+rho*(u-u2);
		if mod(iter,25)==0
            primalcost = sum(sum(sum((x-y).^2)))/2+lambda*nucnorm(opD(x));
            dualcost = cy-sum(sum(sum((y-opDadj(u)).^2)))/2;
				% best value of dualcost computed so far:
			primalcostlowerbound = max(primalcostlowerbound,dualcost);
				% The gap between primalcost and primalcostlowerbound is even better
				% than between primalcost and dualcost to monitor convergence. 
			fprintf('nb iter:%4d  %f  %f  %e\n',iter,primalcost,...
				primalcostlowerbound,primalcost-primalcostlowerbound);
			figure(3);
			imshow(x);
		end
	end
end

% SLOW!
% function val = nucnorm(y)
%     M = size(y,1);
%     N = size(y,2);
%     val = zeros(M,N);
%     for m = 1:M
%         for n = 1:N
%             P_mn = squeeze(y(m,n,:,:));   % 3-by-2
%             [~,S,~] = svd(P_mn,'econ');
%             val(m,n) = sum(diag(S));
%         end
%     end
%     val = sum(val,'all');
% end

function val = nucnorm(y)
    a1 = y(:,:,:,1);
    a2 = y(:,:,:,2);
    a1_norm = sqrt(sum(a1.^2,3));
    a2_norm = sqrt(sum(a2.^2,3));
    theta = atan2(sqrt(sum((cross(a1,a2,3)).^2,3)),dot(a1,a2,3));       % angle between a1 and a2
    val = sqrt(a1_norm.^2 + a2_norm.^2 + 2*sin(theta).*(a1_norm.*a2_norm));
    val = sum(val,'all');
end

% function u = prox_sigma_g_conj(y, lambda)
% 	[M,N,~] = size(y);
%     u = zeros(size(y));
%     for m = 1:M
%         for n = 1:N
%             P_mn = squeeze(y(m,n,:,:));   % 3-by-2
%             [U,S,V] = svd(P_mn,'econ');
%             S_lambda = diag(min(diag(S),lambda));
%             P_mn = U*S_lambda*V';
%             u(m,n,:,:) = reshape(P_mn,1,1,size(P_mn,1),size(P_mn,2));
%         end
%     end
% end

function u = prox_sigma_g_conj(y, lambda)
	s = diff(sum(y.^2,3),1,4);                                  % MxN
	theta = atan2(2*dot(y(:,:,:,1),y(:,:,:,2),3),-s)/2;         % MxN
	c = cos(theta);                                             % MxN
	s = sin(theta);                                             % MxN
    x = cat(4,bsxfun(@times,y(:,:,:,1),c)+bsxfun(@times,y(:,:,:,2),s),...
		bsxfun(@times,y(:,:,:,2),c)-bsxfun(@times,y(:,:,:,1),s));               % MxNx3x2 (U*S)
    S = sqrt(sum(x.^2,3));                      % MxNx1x2 (the norm of the two columns are the two sigular values)
    tmp = min(S, lambda);                       % MxNx1x2 
    tmp = bsxfun(@times, x./(S+eps), tmp);      % MxNx3x2 (each column's norm gets rescaled)
    u = cat(4,bsxfun(@times,tmp(:,:,:,1),c)-bsxfun(@times,tmp(:,:,:,2),s),...
		bsxfun(@times,tmp(:,:,:,2),c)+bsxfun(@times,tmp(:,:,:,1),s));           % MxNx3x2
end
