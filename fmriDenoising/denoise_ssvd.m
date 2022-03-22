function [im_out,rank_out,noise_out] = denoise_ssvd(nim, b, step, wantGaussWeighting,wantRankWeighting,k0)
%
% 'denoise_ssvd' is a wraper of RMT-based denoising algorithm to denoise
% 3D/4D images with the last dimention representing the repetitive
% measurement. Additionally multiple estimates derived from neighboring patches 
% are aggregated for a voxel. The code is modified based on glhosvd.m provided
% by Xinyuan Zhang.
%
% Usage: [im_out, rank_out] = denoise_ssvd (nim, b, step, wantGaussWeighting,wantRankWeighting,k0)
%
% Returns
% -------
% im_out: [x y z M]
%
% rank_out: [x y z], ranks for individual patches after singular value
% shrinkage.
%
% Expects
% -------
% nim: input image [x y z M]
%
% b:   (optional)  window size, defaults to 5 ie 5x5x5 kernel
%
% step: step length between neighboring patches. defaults to 1.
%
% wantGaussWeighting: how much Gaussian weighting is used on the patch.
% Default: [1 400]. The first number is whether to
% use Gaussian weighting and the second number is width factor, inversely
% proportional to the width of the window. 400 means delta weighting.
%
% wantRankWeighting: Whether to use rank weighting, default is 0.
%
% k0: the moments for the generalized circle law to dertermine the rank
% used in ssvd.m for 'ssvd' method.
%
% See also: denoise_mppca denoise_optim_SVHT estimate_noise
%
%
% Copyright (C) 2019 CMRR at UMN
% Author: Xiaoping Wu <xpwu@cmrr.umn.edu> and Wei Zhu <zhuwei@umn.edu>
% Created: Tue Sep  3 15:02:51 2019
%

if ~exist('wantGaussWeighting','var') || isempty(wantGaussWeighting)
  wantGaussWeighting = 1;
end
if ~exist('wantRankWeighting','var') || isempty(wantRankWeighting)
  wantRankWeighting = 0;
end
if ~exist('k0','var') || isempty(k0)
  k0 = 2;
end
if length(wantGaussWeighting)==1
  wantGaussWeighting = [wantGaussWeighting 4];
end


time0         =   clock;
if nargin<2
    b             =   5; % block size
end
if isscalar(b)
    b             =   repmat(b,1,ndims(data)-1); % block size
end
bx = b(1);
by = b(2);
bz = b(3);
if nargin<3
    step          =   1; % step length
end

fprintf('--------start denoising--------\n');
%%%The local denoising stage
[sx,sy,sz,M] = size(nim);
Ys = zeros( size(nim) );
W = zeros( size(nim) );
R = zeros(sx,sy,sz); % rank
S = zeros(sx,sy,sz); % noise standard deviation

%
indices_x= [1:step:sx-bx sx-bx+1];
indices_y= [1:step:sy-by sy-by+1];
indices_z= [1:step:sz-bz sz-bz+1];

% segment the data for parfor
disp('-> segment data...')
data0= zeros(bx,sy,sz,M,length(indices_x));
for i  =  1:length(indices_x) %[1:step:size(nim,1)-b size(nim,1)-b+1]
    data0(:,:,:,:,i)= nim(indices_x(i): indices_x(i)+bx-1, :, :, :);
end

% denoise
disp('-> denoise...')
Ys0= zeros(bx,sy,sz,M,length(indices_x));
W0= Ys0;
R0= zeros(sy,sz,length(indices_x));
S0= R0;

parfor  i  =  1:length(indices_x) %[1:step:size(nim,1)-b size(nim,1)-b+1]
    fprintf('--- denoising: i=%i (%i total) --- \n',i, length(indices_x))
    
    iB1= data0(:, :, :, :,i);
    iYs = zeros(bx,sy,sz,M);
    iW = iYs;
    iR = zeros(sy,sz);
    iS = zeros(sy,sz);
    
    for j = indices_y % [1:step:size(nim,2)-b size(nim,2)-b+1]
        
        %fprintf('--- denoising: i=%i (%i total), j=%i (%i total) --- \n',i, length(indices_x), j, length(indices_y))
        
        for k = indices_z % [1:step:size(nim,3)-b size(nim,3)-b+1]
            
            B1=iB1(:, j:j+by-1, k:k+bz-1, :);
            
            [Ysp, Wp, Rp, Sigmap]   =   Low_rank_SSC(double(B1), wantGaussWeighting, wantRankWeighting,k0);
            
            iYs(:,j:j+by-1,k:k+bz-1,:)=iYs(:,j:j+by-1,k:k+bz-1,:)+ Ysp;
            iW(:,j:j+by-1,k:k+bz-1,:)=iW(:,j:j+by-1,k:k+bz-1,:)+ Wp;
            iR(j,k)= Rp;
            iS(j,k)= Sigmap;
            
        end
    end
    
    Ys0(:,:,:,:,i)= iYs;
    W0(:,:,:,:,i)= iW;
    R0(:,:,i)= iR;
    S0(:,:,i)= iS;
    
end

% aggregate data
disp('-> aggregate segmented results...')
for i  =  1:length(indices_x) %[1:step:size(nim,1)-b size(nim,1)-b+1]
    
    Ys(indices_x(i):indices_x(i)+bx-1, :, :, :)= Ys(indices_x(i):indices_x(i)+bx-1, :, :, :)+ Ys0(:,:,:,:,i);
    W(indices_x(i):indices_x(i)+bx-1, :, :, :)= W(indices_x(i):indices_x(i)+bx-1, :, :, :)+ W0(:,:,:,:,i);
    R(indices_x(i),:,:)= R0(:,:,i);
    S(indices_x(i),:,:)= S0(:,:,i);
    
end

%
im_out  =  Ys./W;
rank_out = R;
noise_out = S;

%Sigma = SIGs./W1;
fprintf('Total elapsed time = %f min\n\n', (etime(clock,time0)/60) );

end

function  [X, W, R, sigma]   =   Low_rank_SSC(Y, wantGaussWeighting,wantRankWeighting,k0)
if ~exist('k0','var') || isempty(k0)
  k0 = 2;
end

[X, R, sigma]= ssvd(Y,'svs1','ssvd',k0);

if wantRankWeighting
    wei = 1/(1+R);
else
    wei = 1;
end

if wantGaussWeighting(1)
    % Weighting changes as ws changes
    wei = wei.*repmat(windowN(@gausswin,sizefull(X,3),wantGaussWeighting(2)),1,1,1,size(X,4));
    
    % weighting fixed for all ws, near to the delta function
    %wei =   wei.*repmat(padarray(windowN(@gausswin,[3 3 3],4),([sx sy sz]-[3 3 3])/2),1,1,1,M);
end

W   =   wei.*ones( size(X) );
X   =   X.*wei;
end






