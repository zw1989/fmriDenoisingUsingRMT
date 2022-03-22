function [imgDenoised,rank,sigma,sigmaVST] = fmriDenoise_SSVD2(epi_noisy, ks, ws, VST_ABC, wantGaussWeighting,wantRankWeighting,wantVST,k0)

%
% 'fmriDenoise_SSVD2' is a wraper of RMT-based denoising algorithm with (or
% without) VST to denoise magnitude (or complex) 3D/4D images with the last
% dimention representing the repetitive measurement. 
%
% Usage: [imgDenoised,rank,sigma,sigmaVST] = fmriDenoise_SSVD2 (epi_noisy, ks, ws, VST_ABC, wantGaussWeighting,wantRankWeighting,wantVST,k0)
%
% Returns
% -------
% imgDenoised: [x y z M]
% rank: [x y z], ranks for individual patches after singular value
%   shrinkage.
% sigma: [x y z], estimated noise
% sigmaVST: [x,y,z], estimated noise from VST
%
% Expects
% -------
% epi_noisy: noisy epi images
% ks: kernel size used in VST operation
% ws: window size, defaults is 3x3x3 
% VST_ABC: which variance stabilizer in VST should be used, 'A' or 'B' (default)
% wantGaussWeighting: how much Gaussian weighting is used on the patch.
%   Default: [1 400]. The first number is whether to
%   use Gaussian weighting and the second number is width factor, inversely
%   proportional to the width of the window. 400 means delta weighting.
% wantRankWeighting: Whether to use rank weighting, default is 0.
%   k0: the moments for the generalized circle law to dertermine the rank
%   used in ssvd.m for 'ssvd' method.
% wantVST: whether use VST to do Rician correction
%
% See also: denoise_ssvd.m ssvd.m
%
%
% Copyright (C) 2019 CMRR at UMN
% Author: Xiaoping Wu <xpwu@cmrr.umn.edu> and Wei Zhu <zhuwei@umn.edu>
% Created: Tue Sep  3 15:02:51 2019

imgSize = size(epi_noisy);

% parameter setting
if ~exist('ks','var') || isempty(ks)
    ks = [min(15,imgSize(1)), min(15,imgSize(2)), min(15,imgSize(3))]; % kernel size for noise estimation
end
if ~exist('ws','var') || isempty(ws)
    ws = [3 3 3]; % window size for denosing
end
if ~exist('VST_ABC','var') || isempty(VST_ABC)
    VST_ABC = 'B'; % choose the function for noise estimation; VST-B recommended
end
if ~exist('wantGaussWeighting','var') || isempty(wantGaussWeighting)
  wantGaussWeighting = 1;
end
if ~exist('wantRankWeighting','var') || isempty(wantRankWeighting)
  wantRankWeighting = 0;
end
if ~exist('wantVST','var') || isempty(wantVST)
  wantVST = 1;
end
if ~exist('k0','var') || isempty(k0)
  k0 = 2;
end


% parallel computation
if isempty(gcp)
    mypool= parpool(6);
end

% figure; imshow(makeMontage(epi0(:,:,sliceSelect,1),1,2,'xy'),[]);
% figure; imshow(makeMontage(epi0_noisy(:,:,sliceSelect,1),1,2,'xy'),[]);


%% If VST is wanted
if wantVST
    % Estimate noise, which method to use??
    sigmaVST = estimate_noise_vst3(epi_noisy,ks,'B') ; 
    % figure; imshow(makeMontage(sigma,1,0,'xy'),[]);

    % VST to Transform noise
    imgRaw= perform_riceVST3(epi_noisy, sigmaVST, ks, VST_ABC) ; 
    
else
    imgRaw = epi_noisy;
    sigmaVST = [];
end

%% denoise using SSVD
step = 1;
[imgDenoised,rank,sigma] = denoise_ssvd(imgRaw,ws,step,wantGaussWeighting,wantRankWeighting,k0);

%% If VST is used, do EUI VST
if wantVST
    % EUI VST
    imgDenoised = perform_riceVST_EUI3(imgDenoised, sigmaVST, ks, VST_ABC);    
 
end




