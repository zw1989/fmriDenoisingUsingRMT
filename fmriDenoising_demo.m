clear; close all; clc;

% This is a demo on denoising fMRI data using RMT-based algorithm with or
% wihtout VST.

%% %%%%%%%%%%%%%%%%%%%%%%% Load EPI Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load all data to cells where data size can vary.

epifilenames = {'sampleFmriData.nii'}; % EDIT THIS

fprintf('loading functional data ...\n')
epis = cell(1,nData); %episize = cell(1,nData); epitr = cell(1,nData);
for p=1:nData
    nii = load_untouch_nii(gunziptemp(epifilenames{p}));
    
    if min(nii.img(:)) < intmin('int16') || max(nii.img(:)) > intmax('int16')
        epis{p} = nii.img/ceil((max(nii.img(:))/double(intmax('int16'))));
    else
        epis{p} = nii.img;     
    end

    % Load some parameters
    if p == 1
        episize = nii.hdr.dime.pixdim(2:4); %mm
        epidim = nii.hdr.dime.dim(2:5);
        epitr = 2; % s
    end 

end
fprintf('done (loading EPI data).\n');

%% EPI volume drop and raw data inspection if necessary


%% Image denoising

fprintf('Denoising raw data... ');
ks = [7 7 3];
ws = [7 7 3];

% Denoise with VST-SSVD, delta weighting(alpha=400), no rank weighting, 
[epis_vstssvd,rank_vstssvd,sigma_vstssvd,sigmaVST] = ...
    cellfun(@(x) fmriDenoise_SSVD2(x, ks, ws, 'A', [], [1 400], 0, 1, 1:10),epis, 'UniformOutput',0);

% Denoise using SSVD without VST
% [epis_ssvd,rank_ssvd,sigma_ssvd] = ...
%     cellfun(@(x) fmriDenoise_SSVD2(x, [10 10 3], ws, 'A', [], [1 400], 0, 0, 1:10), epis, 'UniformOutput',0);

% Save data
cellfun(@(x,y) save_nii(make_nii(x, episize),fullfile(y,['epiVstssvdDeltaWeight_ws',num2str(ws(1)),'.nii'])), epis_vstssvd, epiIntermPath, 'UniformOutput',0); 
cellfun(@(x,y) save_nii(make_nii(x, episize),fullfile(y,['epiVstssvdDeltaWeightSigma_ws',num2str(ws(1)),'.nii'])), sigma_vstssvd, epiIntermPath, 'UniformOutput',0); 
cellfun(@(x,y) save_nii(make_nii(x, episize),fullfile(y,['epiVstssvdDeltaWeightRank_ws',num2str(ws(1)),'.nii'])), rank_vstssvd, epiIntermPath, 'UniformOutput',0); 
cellfun(@(x,y) save_nii(make_nii(x, episize),fullfile(y,['epiVstssvdDeltaWeightSigmavst_ws',num2str(ws(1)),'.nii'])), sigmaVST, epiIntermPath, 'UniformOutput',0); 

fprintf('done.\n'); reportmemoryandtime;


%% Motion estimation and correction 



