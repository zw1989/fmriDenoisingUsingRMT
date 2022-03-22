clear; close all; clc

%% Noise-free matrix X simulation and svd, based on a simple fMRI block design
rng(200);

I = 100; % max image intensity
Xp = round(I*phantom('Modified Shepp-Logan',15));
Xp(sum(Xp,2)==0,:)=[]; % Remove zero rows and columns
Xp(:,sum(Xp,1)==0)=[];
% Xp(Xp==0) = [];
m = numel(Xp);
% figure; imagesc(Xp); axis equal; axis off; colorbar
% figurewrite('%da-SheppLogan',2,{0 [1 300]},figuredir);
% mask = Xp~=0;

beta = 0.55; % aspect ratio m/n
n = fix(m/beta); % e.g. fmri repetation number
% n = 50;
X0 = repmat(Xp(:), 1,n);

duration = 7*n/50+1; tr =2;
hrfknob = getcanonicalhrf(duration,tr); 
design = zeros(1,n); design([10*fix(n/50), 30*fix(n/50)]) = 1;
design = conv2(design,hrfknob);  % convolve design matrix with HRF
design = I*design(1:n); % extract desired subset, I*boldPercent is the bold intensity

% figure; plot(design,'linewidth',5); axis equal; axis tight; axis off;
% figurewrite('%da-stim',2,{0 [1 300]},figuredir);

X0(Xp==100,:) = X0(Xp==100,:) + 0.05*repmat(design,sum(Xp(:)==100),1);
X0(Xp==40,:) = X0(Xp==40,:) + 0.06*repmat(circshift(design,1),sum(Xp(:)==40),1);
X0(Xp==30,:) = X0(Xp==30,:) + 0.01*repmat(circshift(design,2),sum(Xp(:)==30),1);
X0(Xp==20,:) = X0(Xp==20,:) + 0.02*repmat(circshift(design,4),sum(Xp(:)==20),1);
%X0(Xp==10,:) = X0(Xp==10,:) + 0.01*repmat(circshift(design,4),sum(Xp(:)==10),1);

% figure; imagesc(X0); axis equal; axis tight; axis off; colorbar
% figurewrite('%da-SheppLoganFlattened',2,{0 [1 300]},figuredir);

% Centering
% X0c = X0-repmat(mean(X0,2),1,n);

% SVD
[U0,S0,V0] = svd(X0/sqrt(n),'econ'); % Normalize the matrix with sqrt(n)
s0 = diag(S0);
r = length(s0(s0>1E-6)); % Marix rank

%% Denoise Gaussian noise
t = 10; % trials
R_ssvd = zeros(1,t); R_mppca = zeros(1,t); R_gavish = R_ssvd; R_stein = R_ssvd; R_cordero = R_ssvd;
sigma_ssvd = zeros(1,t); sigma_mppca=zeros(1,t); sigma_gavish = sigma_ssvd; sigma_cordero = sigma_ssvd; sigma_stein = sigma_ssvd;
X_ssvd = zeros([size(X0),t]); X_mppca = X_ssvd; X_gavish = X_ssvd; X_stein = X_ssvd; X_cordero = X_ssvd;

mse_ssvd = zeros(1,t);
mse_mppca = mse_ssvd;
mse_gavish = mse_ssvd;
mse_cordero = mse_ssvd;
mse_stein = mse_ssvd;

sigma0 = 1; % Take sigma0 so that sigma0 <= s0(3)/beta^0.25 (the critical point)

rng(200);
for it = 1:t
    fprintf('--- Processing: ii=%i (%i total) --- \n',it, t)
    
    Y = X0 + sigma0*randn(m,n); % Noisy matrix Y = X + epsilon. SNR ~ I/sigma0
    %Y = sqrt((X0+sigma0*randn(m,n)).^2+(sigma0*randn(m,n)).^2);
    
    % figure; imshow(Y, [min(X0(:)), max(X0(:))]);
    % figure; imagesc(Y); axis equal
    
    % Centering
    % Y = Y-repmat(mean(Y,2),1,n);
    
    % My method
    [X_ssvd(:,:,it), R_ssvd(it), sigma_ssvd(it)] = ssvd(Y,'svs1','ssvd',1:10);
    mse_ssvd(it) = sum((svd(X0-X_ssvd(:,:,it))).^2);
    
    % mppca
    [X_mppca(:,:,it), R_mppca(it), sigma_mppca(it)] = ssvd(Y,'none','mppca');
    mse_mppca(it) = sum((svd(X0-X_mppca(:,:,it))).^2);
    
    % gavish
    [X_gavish(:,:,it), R_gavish(it), sigma_gavish(it)] = ssvd(Y,'svs1','gavish');
    mse_gavish(it) = sum((svd(X0-X_gavish(:,:,it))).^2);
    
    % cordero
    [X_cordero(:,:,it), R_cordero(it),sigma_cordero(it)] = ssvd(Y,'svs2','cordero',sigma0,fix(1000/m));
    mse_cordero(it) = sum((svd(X0-X_cordero(:,:,it))).^2);
    
    % stein
    [X_stein(:,:,it), R_stein(it),sigma_stein(it)] = ssvd(Y,'none','stein',sigma0*randn(5*m,5*n),20);
    mse_stein(it) = sum((svd(X0-X_stein(:,:,it))).^2);

    
end

%% Plot  
R_all = [R_ssvd', R_mppca',R_cordero',R_stein'];
sigma_all = [sigma_ssvd', sigma_mppca',sigma_cordero',sigma_stein'];
mse_all = [mse_ssvd', mse_mppca',mse_cordero',mse_stein'];

xLabel = {'SVS-proposed','MPPCA','SVS-Cordero','NORDIC'}; 
xLabel = {'','','',''};

figureprep([100 100 1300 400],1); 
ha = tight_subplot(1,3,[.07 .07],[.2 .1],[.05 .05]);

% Plot estimated rank
subplot('Position',ha(1,:)); 
h = iosr.statistics.boxPlot(xLabel,R_all,...
  'symbolMarker','*','symbolColor',[0.85,0.3,0.2],'medianColor',[0.85,0.3,0.2],...
  'lineWidth',2,'lineColor',[0.2 0.3 0.95],'lineStyle','-',...
  'boxcolor',[0.9 0.95 0.9],'boxWidth', 0.5,...
  'showViolin', false, 'violinColor',[0.8,0.8,0.8],...
  'showMean',true, 'meanColor','k', 'meanMarker','.','meanSize',20);
box on; ylabel('Estimated rank'); % xlabel('Methods');
%ylim([1.8, 4.2]);
h.handles.axes.XTickLabelRotation = 10;

hold on; plot(0:6,r*ones(7,1),'k:', 'linewidth',1); hold off% Plot true rank
mean(R_all)

% Plot estimated noise
subplot('Position',ha(2,:)); 
h = iosr.statistics.boxPlot(xLabel,sigma_all,...
  'symbolMarker','*','symbolColor',[0.85,0.3,0.2],'medianColor',[0.85,0.3,0.2],...
  'lineWidth',2,'lineColor',[0.2 0.3 0.95],'lineStyle','-',...
  'boxcolor',[0.9 0.95 0.9],'boxWidth', 0.5,...
  'showViolin', false, 'violinColor',[0.8,0.8,0.8],...
  'showMean',true, 'meanColor','k', 'meanMarker','.','meanSize',20);
box on; ylabel('Estimated noise');  % xlabel('Methods');
% ylim([1.8, 4.2]);
h.handles.axes.XTickLabelRotation = 10;
h.handles.axes.YTick = 0.94:0.02:1.03;

hold on; plot(0:6,sigma0*ones(7,1),'k:', 'linewidth',1); hold off% Plot true rank
mean(sigma_all)

% Plot MSE
subplot('Position',ha(3,:));
h = iosr.statistics.boxPlot(xLabel,mse_all,...
  'symbolMarker','*','symbolColor',[0.85,0.3,0.2],'medianColor',[0.85,0.3,0.2],...
  'lineWidth',2,'lineColor',[0.2 0.3 0.95],'lineStyle','-',...
  'boxcolor',[0.9 0.95 0.9],'boxWidth', 0.5,...
  'showViolin', false, 'violinColor',[0.8,0.8,0],...
  'showMean',true, 'meanColor','k', 'meanMarker','.','meanSize',20);
box on
ylabel('MSE'); % xlabel('Methods');
h.handles.axes.XTickLabelRotation = 10;
h.handles.axes.YTick = 500:500:2500;

h = findobj(gcf); % Get the handles associated with the current figure
allAxes = findall(h, 'Type', 'axes'); set(allAxes, 'Linewidth',1, 'FontWeight','bold', 'FontSize',16, 'Box', 'on',...
            'xgrid','off','ygrid','off','GridLineStyle','--','GridColor',[0.005 0.005,0.005],'XTickLabelRotation',11);
allText = findall(h, 'Type', 'text'); set(allText, 'FontWeight','bold', 'FontSize',24)

% figurewrite('%da-pcaMethodComparison-largerFont-noxlabel',4,{0 [1 300]},figuredir);


% figure; imshow(X_ssvd,[])
% figure; imshow(X_mppca,[])
% figure; imshow(X_cordero,[])
% figure; imshow(X0c,[])


%% Denoise Rician noise
t = 10; % trials
R_ssvd = zeros(1,t); R_mppca = zeros(1,t); R_gavish = R_ssvd; R_stein = R_ssvd; R_cordero = R_ssvd;
sigma_ssvd = zeros(1,t); sigma_mppca=zeros(1,t); sigma_gavish = sigma_ssvd; sigma_cordero = sigma_ssvd; sigma_stein = sigma_ssvd;
X_ssvd = zeros([size(X0),t]); X_mppca = X_ssvd; X_gavish = X_ssvd; X_stein = X_ssvd; X_cordero = X_ssvd;

mse_ssvd = zeros(1,t);
mse_mppca = mse_ssvd;
mse_gavish = mse_ssvd;
mse_cordero = mse_ssvd;
mse_stein = mse_ssvd;

sigma0 = 1; % Take sigma0 so that sigma0 <= s0(3)/beta^0.25 (the critical point)

rng(200);
for it = 1:t
    fprintf('--- Processing: ii=%i (%i total) --- \n',it, t)
    
    %Y = X0 + sigma0*randn(m,n); % Noisy matrix Y = X + epsilon. SNR ~ I/sigma0
    Y = sqrt((X0+sigma0*randn(m,n)).^2+(sigma0*randn(m,n)).^2);
    
    % figure; imshow(Y, [min(X0(:)), max(X0(:))]);
    % figure; imagesc(Y); axis equal
    
    % Centering
    % Y = Y-repmat(mean(Y,2),1,n);
    
    % My method
    [X_ssvd(:,:,it), R_ssvd(it), sigma_ssvd(it)] = ssvd(Y,'svs1','ssvd',1:10);
    mse_ssvd(it) = sum((svd(X0-X_ssvd(:,:,it))).^2);
    
    % mppca
    [X_mppca(:,:,it), R_mppca(it), sigma_mppca(it)] = ssvd(Y,'none','mppca');
    mse_mppca(it) = sum((svd(X0-X_mppca(:,:,it))).^2);
    
    % gavish
    [X_gavish(:,:,it), R_gavish(it), sigma_gavish(it)] = ssvd(Y,'svs1','gavish');
    mse_gavish(it) = sum((svd(X0-X_gavish(:,:,it))).^2);
    
    % cordero
    [X_cordero(:,:,it), R_cordero(it),sigma_cordero(it)] = ssvd(Y,'svs2','cordero',sigma0,fix(1000/m));
    mse_cordero(it) = sum((svd(X0-X_cordero(:,:,it))).^2);
    
    % stein
    [X_stein(:,:,it), R_stein(it),sigma_stein(it)] = ssvd(Y,'none','stein',sigma0*randn(5*m,5*n),20);
    mse_stein(it) = sum((svd(X0-X_stein(:,:,it))).^2);

    
end

%% Plot  
R_all = [R_ssvd', R_mppca',R_cordero',R_stein'];
sigma_all = [sigma_ssvd', sigma_mppca',sigma_cordero',sigma_stein'];
mse_all = [mse_ssvd', mse_mppca',mse_cordero',mse_stein'];

xLabel = {'SVS-proposed','MPPCA','SVS-Cordero','NORDIC'}; 
xLabel = {'','','',''};

figureprep([100 100 1300 400],1); 
ha = tight_subplot(1,3,[.07 .07],[.2 .1],[.05 .05]);

% Plot estimated rank
subplot('Position',ha(1,:)); 
h = iosr.statistics.boxPlot(xLabel,R_all,...
  'symbolMarker','*','symbolColor',[0.85,0.3,0.2],'medianColor',[0.85,0.3,0.2],...
  'lineWidth',2,'lineColor',[0.2 0.3 0.95],'lineStyle','-',...
  'boxcolor',[0.9 0.95 0.9],'boxWidth', 0.5,...
  'showViolin', false, 'violinColor',[0.8,0.8,0.8],...
  'showMean',true, 'meanColor','k', 'meanMarker','.','meanSize',20);
box on; ylabel('Estimated rank'); % xlabel('Methods');
%ylim([1.8, 4.2]);
h.handles.axes.XTickLabelRotation = 10;

hold on; plot(0:6,r*ones(7,1),'k:', 'linewidth',1); hold off% Plot true rank
mean(R_all)

% Plot estimated noise
subplot('Position',ha(2,:)); 
h = iosr.statistics.boxPlot(xLabel,sigma_all,...
  'symbolMarker','*','symbolColor',[0.85,0.3,0.2],'medianColor',[0.85,0.3,0.2],...
  'lineWidth',2,'lineColor',[0.2 0.3 0.95],'lineStyle','-',...
  'boxcolor',[0.9 0.95 0.9],'boxWidth', 0.5,...
  'showViolin', false, 'violinColor',[0.8,0.8,0.8],...
  'showMean',true, 'meanColor','k', 'meanMarker','.','meanSize',20);
box on; ylabel('Estimated noise');  % xlabel('Methods');
% ylim([1.8, 4.2]);
h.handles.axes.XTickLabelRotation = 10;
h.handles.axes.YTick = 0.6:0.05:1.03;

hold on; plot(0:6,sigma0*ones(7,1),'k:', 'linewidth',1); hold off% Plot true rank
mean(sigma_all)

% Plot MSE
subplot('Position',ha(3,:));
h = iosr.statistics.boxPlot(xLabel,mse_all,...
  'symbolMarker','*','symbolColor',[0.85,0.3,0.2],'medianColor',[0.85,0.3,0.2],...
  'lineWidth',2,'lineColor',[0.2 0.3 0.95],'lineStyle','-',...
  'boxcolor',[0.9 0.95 0.9],'boxWidth', 0.5,...
  'showViolin', false, 'violinColor',[0.8,0.8,0],...
  'showMean',true, 'meanColor','k', 'meanMarker','.','meanSize',20);
box on
ylabel('MSE'); % xlabel('Methods');
h.handles.axes.XTickLabelRotation = 10;
h.handles.axes.YTick = 1.2E4:1000:1.6E4;

h = findobj(gcf); % Get the handles associated with the current figure
allAxes = findall(h, 'Type', 'axes'); set(allAxes, 'Linewidth',1, 'FontWeight','bold', 'FontSize',16, 'Box', 'on',...
            'xgrid','off','ygrid','off','GridLineStyle','--','GridColor',[0.005 0.005,0.005],'XTickLabelRotation',11);
allText = findall(h, 'Type', 'text'); set(allText, 'FontWeight','bold', 'FontSize',24)

% figurewrite('%db-pcaMethodComparison-Rician-largerFont-noxlabel',4,{0 [1 300]},figuredir);


%% Denoise Racian noise After VST
t = 10; % trials
R_ssvd = zeros(1,t); R_mppca = zeros(1,t); R_gavish = R_ssvd; R_stein = R_ssvd; R_cordero = R_ssvd;
sigma_ssvd = zeros(1,t); sigma_mppca=zeros(1,t); sigma_gavish = sigma_ssvd; sigma_cordero = sigma_ssvd; sigma_stein = sigma_ssvd;
X_ssvd = zeros([size(X0),t]); X_mppca = X_ssvd; X_gavish = X_ssvd; X_stein = X_ssvd; X_cordero = X_ssvd;

mse_ssvd = zeros(1,t);
mse_mppca = mse_ssvd;
mse_gavish = mse_ssvd;
mse_cordero = mse_ssvd;
mse_stein = mse_ssvd;

sigma0 = 1; % Take sigma0 so that s0(2)/sigma0 <= beta^0.25 (the critical point)

rng(200);
for it = 1:t
    fprintf('--- Processing: ii=%i (%i total) --- \n',it, t)
    
    % Y = X0 + sigma0*randn(m,n); % Noisy matrix Y = X + epsilon. SNR ~ I/sigma0
    Y = sqrt((X0+sigma0*randn(m,n)).^2+(sigma0*randn(m,n)).^2);
    % sigma_tmp = riceVST_sigmaEst(Y,0,1,'B');
    sigma_tmp = sigma0;
    %[~,~, sigma_tmp] = ssvd(Y,'none','gavish');
    Y = riceVST(Y,sigma_tmp,'A'); % VST to Transform noise
    % figure; imshow(Y, [min(X0(:)), max(X0(:))]);
    % figure; imagesc(Y); axis equal
    
    % Centering
    % Y = Y-repmat(mean(Y,2),1,n);
    
    % My method
    [X_ssvd(:,:,it), R_ssvd(it), sigma_ssvd(it)] = ssvd(Y,'svs1','ssvd',1:10);
    X_ssvd(:,:,it) = riceVST_EUI(X_ssvd(:,:,it), sigma_tmp, 'A');
    mse_ssvd(it) = sum((svd(X0-X_ssvd(:,:,it))).^2);
    
    % mppca
    [X_mppca(:,:,it), R_mppca(it), sigma_mppca(it)] = ssvd(Y,'none','mppca');
    X_mppca(:,:,it) = riceVST_EUI(X_mppca(:,:,it), sigma_tmp, 'A');
    mse_mppca(it) = sum((svd(X0-X_mppca(:,:,it))).^2);
    
    % gavish
    [X_gavish(:,:,it), R_gavish(it), sigma_gavish(it)] = ssvd(Y,'svs1','gavish');
    X_gavish(:,:,it) = riceVST_EUI(X_gavish(:,:,it), sigma_tmp, 'A');
    mse_gavish(it) = sum((svd(X0-X_gavish(:,:,it))).^2);
    
    % cordero
    [X_cordero(:,:,it), R_cordero(it),sigma_cordero(it)] = ssvd(Y,'svs2','cordero',sigma_tmp,fix(1000/m));
    X_cordero(:,:,it) = riceVST_EUI(X_cordero(:,:,it), sigma_tmp, 'A');
    mse_cordero(it) = sum((svd(X0-X_cordero(:,:,it))).^2);
    
    % stein
    [X_stein(:,:,it), R_stein(it),sigma_stein(it)] = ssvd(Y,'none','stein',sigma_tmp*randn(5*m,5*n),20);
    X_stein(:,:,it) = riceVST_EUI(X_stein(:,:,it), sigma_tmp, 'A');
    mse_stein(it) = sum((svd(X0-X_stein(:,:,it))).^2);

    
end

%% Plot  
R_all = [R_ssvd', R_mppca',R_cordero',R_stein'];
sigma_all = [sigma_ssvd', sigma_mppca',sigma_cordero',sigma_stein'];
mse_all = [mse_ssvd', mse_mppca',mse_cordero',mse_stein'];

xLabel = {'SVS-proposed','MPPCA','SVS-Cordero','NORDIC'}; 
xLabel = {'','','',''};

figureprep([100 100 1300 400],1); 
ha = tight_subplot(1,3,[.07 .07],[.2 .1],[.05 .05]);

% Plot estimated rank
subplot('Position',ha(1,:)); 
h1 = iosr.statistics.boxPlot(xLabel,R_all,...
  'symbolMarker','*','symbolColor',[0.85,0.3,0.2],'medianColor',[0.85,0.3,0.2],...
  'lineWidth',2,'lineColor',[0.2 0.3 0.95],'lineStyle','-',...
  'boxcolor',[0.9 0.95 0.9],'boxWidth', 0.5,...
  'showViolin', false, 'violinColor',[0.8,0.8,0.8],...
  'showMean',true, 'meanColor','k', 'meanMarker','.','meanSize',20);
box on; ylabel('Estimated rank'); % xlabel('Methods');
%ylim([1.8, 4.2]);
h1.handles.axes.XTickLabelRotation = 10;

hold on; plot(0:6,r*ones(7,1),'k:', 'linewidth',1); hold off% Plot true rank
mean(R_all)

% Plot estimated noise
subplot('Position',ha(2,:)); 
h2 = iosr.statistics.boxPlot(xLabel,sigma_all,...
  'symbolMarker','*','symbolColor',[0.85,0.3,0.2],'medianColor',[0.85,0.3,0.2],...
  'lineWidth',2,'lineColor',[0.2 0.3 0.95],'lineStyle','-',...
  'boxcolor',[0.9 0.95 0.9],'boxWidth', 0.5,...
  'showViolin', false, 'violinColor',[0.8,0.8,0.8],...
  'showMean',true, 'meanColor','k', 'meanMarker','.','meanSize',20);
box on; ylabel('Estimated noise');  % xlabel('Methods');
% ylim([1.8, 4.2]);
h2.handles.axes.XTickLabelRotation = 10;
h2.handles.axes.YTick = 0.6:0.05:1.05;

hold on; plot(0:6,sigma0*ones(7,1),'k:', 'linewidth',1); hold off% Plot true rank
mean(sigma_all)

% Plot MSE
subplot('Position',ha(3,:));
h3 = iosr.statistics.boxPlot(xLabel,mse_all,...
  'symbolMarker','*','symbolColor',[0.85,0.3,0.2],'medianColor',[0.85,0.3,0.2],...
  'lineWidth',2,'lineColor',[0.2 0.3 0.95],'lineStyle','-',...
  'boxcolor',[0.9 0.95 0.9],'boxWidth', 0.5,...
  'showViolin', false, 'violinColor',[0.8,0.8,0],...
  'showMean',true, 'meanColor','k', 'meanMarker','.','meanSize',20);
box on
ylabel('MSE'); % xlabel('Methods');
h3.handles.axes.XTickLabelRotation = 10;
h3.handles.axes.YTick = 1000:500:4000;

h = findobj(gcf); % Get the handles associated with the current figure
allAxes = findall(h, 'Type', 'axes'); set(allAxes, 'Linewidth',1, 'FontWeight','bold', 'FontSize',16, 'Box', 'on',...
            'xgrid','off','ygrid','off','GridLineStyle','--','GridColor',[0.005 0.005,0.005],'XTickLabelRotation',11);       
allText = findall(h, 'Type', 'text'); set(allText, 'FontWeight','bold', 'FontSize',24)


% figurewrite('%db-pcaMethodComparison-vst-largerFont-noxlabel',4,{0 [1 300]},figuredir);



