function fz_data = perform_riceVST3 (data, sigma, b, VST_ABC)
% PERFORM_RICEVST3 perform riceVST method (ie, variance stabilization
% transform) for spatially varying noise levels. 
% The multiple estimates of a voxel obtained from the sliding window 
% is aggregated to derive an average estimate for that voxel.
% The average noise level within the sliding window is used for VST.
% this is basically a modification of perform_riceVST.m to enable parfor.
%
% Usage: fz_data = perform_riceVST3 (data, sigma, b, VST_ABC)
%
% Returns
% -------
% fz_data: [x, y, z, M] transformed data
%
% Expects
% -------
% data: [x, y, z, M] data matrix
% 
% sigma: [x, y, z] noise map
% 
% b: (optional)  window size, defaults to 5 ie 5x5x5 kernel
% 
% VST_ABC: name or filename of variance-stabilizing transform to be used (default='B')
%
%
% See also: estimate_noise_vst3.m perform_riceVST2.m riceVST.m perform_riceVST_EUI.m riceVST_EUI.m 
%
%
% Copyright (C) 2019 CMRR at UMN
% Author: Xiaoping Wu <xpwu@cmrr.umn.edu> 
% Created: Tue Sep  3 14:34:04 2019
%

if isa(data,'integer')
    data = single(data);
end

if ~exist('VST_ABC', 'var') || isempty(VST_ABC)
    VST_ABC= 'B';
end

time0         =   clock;
if nargin<2
    b             =   5; % block size
end
if isscalar(b)
    b             =   repmat(b,1,ndims(data)-1); % block size
end
%if nargin<3
step          =   1; % step length
%end

fprintf('--------start VST --------\n');
%%%The local mppca denoising stage
Ys            =   zeros( size(data) );
W             =   Ys;
bx = b(1);
by = b(2);
bz = b(3);

len_i= length([1:step:size(data,1)-bx size(data,1)-bx+1]);
len_j= length([1:step:size(data,2)-by size(data,2)-by+1]);

% segment the data for parfor
disp('-> segment data...')
[~,sy,sz,M] = size(data);
data0= zeros(bx,sy,sz,M,len_i);
sigma0= zeros(bx,sy,sz,len_i);
for i  =  [1:step:size(data,1)-bx size(data,1)-bx+1]
    data0(:,:,:,:,i)= data(i:i+bx-1, :, :, :);
    sigma0(:,:,:,i)= sigma(i:i+bx-1, :, :);
end

% VST
disp('-> run VST...')
Ys0= zeros(bx,sy,sz,M,len_i);
W0= Ys0;
parfor  i  =  [1:step:size(data,1)-bx size(data,1)-bx+1]
     fprintf('--- i=%i (%i total) --- \n',i, len_i)
    
    iB1= data0(:, :, :, :,i);
    iSigma= sigma0(:,:,:,i);
    iYs = zeros(bx,sy,sz,M);
    iW= iYs;
    
    for j = [1:step:size(data,2)-by size(data,2)-by+1]
        
        % fprintf('--- VST: i=%i (%i total), j=%i (%i total) --- \n',i, len_i, j, len_j)
        
        
        for k = [1:step:size(data,3)-bz size(data,3)-bz+1]
            
            %B1=data(i:i+b-1, j:j+b-1, k:k+b-1, :);
            B1=iB1(:, j:j+by-1, k:k+bz-1, :);
            
            %Sig= sigma(i:i+b-1, j:j+b-1, k:k+b-1);
            Sig= iSigma(:, j:j+by-1, k:k+bz-1);
            
            rimavst = riceVST(B1,mean(Sig(:)),VST_ABC);  
            
%             Ys(i:i+b-1,j:j+b-1,k:k+b-1,:)=Ys(i:i+b-1,j:j+b-1,k:k+b-1,:)+ rimavst;
%             W(i:i+b-1,j:j+b-1,k:k+b-1,:)=W(i:i+b-1,j:j+b-1,k:k+b-1,:)+ 1;
            
            iYs(:,j:j+by-1,k:k+bz-1,:)=iYs(:,j:j+by-1,k:k+bz-1,:)+ rimavst;
            iW(:,j:j+by-1,k:k+bz-1,:)=iW(:,j:j+by-1,k:k+bz-1,:)+ 1;
            
        end
    end
    
    Ys0(:,:,:,:,i)= iYs;
    W0(:,:,:,:,i)= iW;
    
end
% aggregate data
disp('-> aggregate segmented results...')
for i  =  [1:step:size(data,1)-bx size(data,1)-bx+1]
    
    Ys(i:i+bx-1, :, :, :)= Ys(i:i+bx-1, :, :, :)+ Ys0(:,:,:,:,i);
    W(i:i+bx-1, :, :, :)= W(i:i+bx-1, :, :, :)+ W0(:,:,:,:,i);
    
end

%
fz_data  =  Ys./W;
fprintf('Total elapsed time = %f min\n\n', (etime(clock,time0)/60) );

end


