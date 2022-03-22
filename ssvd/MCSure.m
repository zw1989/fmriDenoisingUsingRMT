%% MC SURE algorithm to estimate MSE without ground truth
%% Reference Ramani, S., Blu, T., Unser, M., 2008. Monte-Carlo Sure: A black-box optimization of regularization parameters for general denoising algorithms. IEEE Transactions on Image Processing 17, 1540-1554.
function SURE = MCSure(y,f,epsilon,sigma,bprime)
%% Input:
% y: noisy image (Nx*Ny)
% f: denoising function, should be pre-defined as: f = @denoise_fun
% epsilon: desired divergence (recomended: 0.01)
% sigma: std of the zero-mean white Gaussian noise on y compared to noise-free image
% bprime: a zero-mean i.i.d. random vector (Nx*Ny) with unit variance and bounded higher order moments
%% Output:
% SURE: the estimated MSE
%% example:
% Img0 = phantom(64,64); % noise-free img
% sigma = 0.05; 
% Img = Img0 + sigma * randn(size(Img0)); % noisy image
% h = fspecial('gaussian',[3,3],0.4); 
% f = @(x)imfilter(x,h);
% epsilon = 0.01;
% bprime = randn(size(Img0));
% SURE = MCSure(Img,f,epsilon,sigma,bprime);
% MSE = sum((Img0(:)-Img(:)).^2)/size(Img,1);
%%
fy = f(y); % denoised image
z = y + epsilon .* bprime;
fz = f(z);
[M,N] = size(y);
div = 1/epsilon * bprime(:)' * (fz(:) - fy(:));
SURE = 1/N*norm(y(:)-fy(:))^2 - M*sigma^2 + 2*sigma^2/N*div;
end
