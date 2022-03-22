% This function creates a N-dimensional window for a sample image, it takes
% the dimensions of the window and applies the 1D window function. It is
% simply an extenstion of the function window in matlab.
% 
% Usage: w = windowN(wfunc,winSize)
%
% Inputs: 
% <winSize> is an integer array with the size of a window you created, the
%   length of which is the window dimension. The window can be any dimension.
% <wfunc> is the functional handle the same as the input for window
%   function. For example, if you want to construct a hamming window, set
%   wfunc = @hamming
%
% Output:
% <w> is the constructed multi-dimensional window
% 
% Example1:
% w = windowN(@hamming, [30,40]);
% figure; imagesc(w)
%
% Example2:
% w = windowN({@hamming,@gausswin}, [30,40]);
% figure; imagesc(w)
%
% By Wei Zhu, zhuwei@umn.edu


function w = windowN(wfunc,winSize, varargin)

% Dimension of the window
ndim = length(winSize);

% Check input
if ~iscell(wfunc)
    wfunc = {wfunc};
end
if length(wfunc)==1
    wfunc = repmat(wfunc,1,ndim); 
end

% Construct windows along each dimension, the size may vary
wtmp = cellfun(@(x,y) window(x,y,varargin{:}), wfunc, num2cell(winSize),'UniformOutput',0);

% Create N-dimensional meshgrid
masktmp = cell(1,ndim);
[masktmp{:}] = ndgrid(wtmp{:}); 

% Final window
w = 1;
for idim = 1:ndim
    w = w.*masktmp{idim};
end

end