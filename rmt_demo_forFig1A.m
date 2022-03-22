clear; close all; clc

%% Noise-free matrix X 
% Note here X is constructed based on a simple fMRI block design
% experiment. In real test, the simplest X construction is a diagonal
% matrix. The number of the non-zero diagonal elements is the rank.

rng(200);

% Simulation matrix m x n, assume m <= n
beta = 0.55; % aspect ratio m/n
n = 300; % e.g. fmri repetation number
m = round(n*beta); % e.g. voxel number in a 3D patch

% Construct noise-free matrix: anatomical structure and BOLD response
X = zeros(m,n);
mp = randi(m); % Voxel number in the anatomical structure
mIdx = randperm(m,mp); % Voxel index
X(mIdx,:) = X(mIdx,:) + 1*repmat(0.8+0.4*rand(mp,1),1,n); % Anatomical structure, in reality the intensity should be high, for demonstration, we take it to be 1

% Construct design matrix
duration = 7*n/50+1; tr =2;
hrfknob = getcanonicalhrf(duration,tr); 
design = zeros(1,n); design([10*n/50, 30*n/50]) = 1;
design = conv2(design,hrfknob);  % convolve design matrix with HRF
design = 4*design(1:n); % extract desired subset, 4 is the signal intensity
% figure; plot(design)

% Construct BOLD signal
np = randi(mp); % How many voxels have bold response?
nIdx = randperm(mp,np);
X(mIdx(nIdx),:) = X(mIdx(nIdx),:) + repmat(design,np,1).*repmat(abs(randn(np,1)),1,n);

% Introduce more ranks
% r0 = 3; % rank of the true matrix
% tmp = beta^0.25 + 10;
% X(1,1) = X(1,1) + tmp*sqrt(n);

figure; imshow(X,[min(X(:)), max(X(:))/2]); %axis on

% SVD
[U0,S0,V0] = svd(X/sqrt(n),'econ'); % Normalize the matrix with sqrt(n)
s0 = diag(S0);
r = length(s0(s0>1E-6)); % Marix rank

% figure; imshow(U0,[-1,1])
% figure; imshow(S0,[])
% figure; imshow(V0',[-1,1])
% figure; scatter(U0(:,1),U0(:,2)); xlabel('PC 1'); ylabel('PC 2'); grid on
% figure; plot(s0)

% Plot singular value histogram
figure; hold on
histogram(s0, linspace(min(s0),max(s0),n),'normalization','pdf'); %title('Singular value density distribution');
ylim([0,0.4]); ylabel('Density'); xlabel('Singular value (SV)'); 
L = legend('Ground truth SV density'); legend('boxoff');
set(L,'Position',[0.4 0.78 0.02 0.02])


%% Noisy matrix Y = X + epsilon. Phenomena2: Eigenvalue bias, largest SV distribution

t = 1000; % trials
s_all = []; % singular value samples
u1_all = []; v1_all = [];
u2_all = []; v2_all = [];
sigma0 = 8; % Take sigma0 so that sigma0 <= s0(2)/beta^0.25 (the critical point)

% Perform experiment
for ii = 1:t
    Y = X + sigma0*randn(m,n);
    %Y = sqrt((1*X+sigma0*randn(m,n)).^2+(sigma0*randn(m,n)).^2);
    %Y = riceVST(Y,sigma0,'B'); % VST to Transform noise

    [U,S,V] = svd(Y/sqrt(n),'econ');
    s = diag(S);
    s_all = [s_all,s]; 
    u1_all = [u1_all, U(:,1)];
    u2_all = [u2_all, U(:,2)];
    v1_all = [v1_all, V(:,1)];
    v2_all = [v2_all, V(:,2)];
end

%% Plot singular value histogram
figureprep([100 100 1600 900],1); 
h = histogram(s_all, linspace(min(s_all(:)),max(s_all(:)),n),'normalization','pdf'); % Sample SV density
histogram([sigma0*ones((m-r)*t,1);repmat(s0(s0>1E-6),t,1)], linspace(min(s_all(:)),max(s_all(:)),n),'normalization','pdf'); % Population SV density
ylabel('Density'); xlabel('Singular value (SV)'); ylim([0,0.12]); %xlim([0,22])

% Theory--PDF for singular value distribution
% Note the integral of the second order moment of the distribution is variance
beta0 = (m-r)/(n-r);
beta1s = (1-sqrt(beta0))*sigma0; beta2s = (1+sqrt(beta0))*sigma0; % Bulk edges for singluar value, theory
ss = linspace(beta1s, beta2s);
qc = @(x,beta, sigma) sqrt((x.^2-(sigma*(1-sqrt(beta))).^2) .* ((sigma*(1+sqrt(beta))).^2-x.^2))./(pi*x*beta*sigma^2); % Quarter-circle distribution
plot(ss, qc(ss, beta0, sigma0), 'linewidth',2); ylabel('Probability'); xlabel('Singular value (SV)'); 

% Theory -- transition point
transitionPoint = beta^0.25*sigma0;
plot(transitionPoint,0,'r.', 'MarkerSize', 40);

% Theory -- Signal SV distribution (Gaussian)
snr = s0(1)/sigma0;
rho1 = sqrt( (1+snr^2) * (beta+snr^2) / snr^2 );
a = integral(@(x) qc(x,beta,1)./(rho1.^2-x.^2).^2, (1-sqrt(beta)), (1+sqrt(beta)))./...
    (integral(@(x) qc(x,beta,1)./(rho1.^2-x.^2), (1-sqrt(beta)), (1+sqrt(beta)))).^2;
b = (beta*integral(@(x) qc(x,beta,1)./(rho1.^2-x.^2).^2, (1-sqrt(beta)), (1+sqrt(beta))) + (1-beta)/rho1^4)./...
    (beta*integral(@(x) qc(x,beta,1)./(rho1.^2-x.^2), (1-sqrt(beta)), (1+sqrt(beta))) + (1-beta)/rho1^2).^2;
c = 2*integral(@(x) qc(x,beta,1).*x.^2./(rho1.^2-x.^2).^2, (1-sqrt(beta)), (1+sqrt(beta)))./...
    integral(@(x) qc(x,beta,1).*rho1./(rho1.^2-x.^2), (1-sqrt(beta)), (1+sqrt(beta)))./...
    (beta*integral(@(x) qc(x,beta,1).*rho1./(rho1.^2-x.^2), (1-sqrt(beta)), (1+sqrt(beta))) + (1-beta)/rho1);

f2 = a+b+c;
tau1 = sqrt(f2/2) * n^(-1/2);

x = -3:0.01:3;
s1_distribution = 1/sqrt(2*pi).*exp(-x.^2/2); 
s1_distribution = s1_distribution/(tau1*sigma0)/m; % Normalize. The idea is to keep the integral of the streched distribution to be still 1.
plot((x*tau1 + rho1)*sigma0, s1_distribution, 'color',[0.2 0.8 0], 'linewidth',2)

% Theory -- noise SV or small signal SV distribution (Tracy-Widom distribution)
tau2 = (1+sqrt(beta))^(4/3) * beta^(-1/6) * n^(-2/3);
rho2 = (1+sqrt(beta))^2;
x = -6:0.001:5;
s2_distribution = tracywidomDistribution(x,1);
s2_distribution = s2_distribution/(1/2*tau2/sqrt(rho2)*sigma0)/m; % Normalize. Here the x value will be changed non-linearly. Tylor expansion on the sqrt(x+rho2/tau2) is used.
plot(sqrt(flip(x)*tau2 + rho2)*sigma0, s2_distribution,'color',[0.2 0 1],'linewidth',2)


L = legend('Sample SV density','Population SV density','Generalized quarter circle law',...
           'Transition point','Gaussian distribution', 'Tracy-Widom distribution'); legend('boxoff');
% set(L,'Position',[0.6 0.85 0.02 0.02]); 
set(L,'Position',[0.7 0.75 0.02 0.02]); 

hold off

h = findobj(gcf); % Get the handles associated with the current figure
allLines = findall(h, 'Type', 'line'); set(allLines, 'Linewidth',4)
allAxes = findall(h, 'Type', 'axes'); set(allAxes, 'Linewidth',2, 'FontWeight','bold', 'FontSize',32, 'Box', 'off',...
            'xgrid','off','ygrid','off','GridLineStyle','--','GridColor',[0.005 0.005,0.005]); 
allText = findall(h, 'Type', 'text'); set(allText, 'FontWeight','bold', 'FontSize',32)

% figurewrite('%da-rmt-largerFont',1,{0 [1 300]},figuredir);

%% Plot singular vectors
uu1 = (u1_all'*U0(:,1)).^2; 
uu10 = 1-beta*(1+(s0(1)/sigma0).^2)./(s0(1)/sigma0).^2./((s0(1)/sigma0).^2+beta);
theta1 = acos(sqrt(uu1))/pi*180;
theta10 = acos(sqrt(uu10))/pi*180;
figure;hold on
histogram(theta1); plot(theta10,0, 'r.', 'MarkerSize', 20); hold off

uu2 = u2_all'*U0(:,2);
uu20 = 0;
% uu20 = 1-beta*(1+(s0(2)/sigma0).^2)./(s0(2)/sigma0).^2./((s0(2)/sigma0).^2+beta);
theta2 = acos(uu2)/pi*180;
% theta2(s_all(2,:)<=beta2s) = [];
theta20 = acos(sqrt(uu20))/pi*180;
figure;hold on
histogram(theta2); plot(theta20,0, 'r.', 'MarkerSize', 20); hold off

vv1 = (v1_all'*V0(:,1)).^2; 
vv10 = 1-(beta+(s0(1)/sigma0).^2)./(s0(1)/sigma0).^2./((s0(1)/sigma0).^2+1);
alpha1 = acos(sqrt(vv1))/pi*180;
alpha10 = acos(sqrt(vv10))/pi*180;
figure;hold on
histogram(alpha1); plot(alpha10,0, 'r.', 'MarkerSize', 20); hold off

vv2 = v2_all'*V0(:,2);
vv20 = zeros(t,1);
alpha2 = acos(vv2)/pi*180;
alpha20 = acos(sqrt(vv20))/pi*180;
figure;hold on
histogram(alpha2); plot(alpha20,0, 'r.', 'MarkerSize', 20); hold off


% Generate the x,y,z data for a sphere
radius0 = 1*ones(t,t); % Radius is 1
[th0, phi0] = meshgrid(linspace(0,2*pi,t), linspace(-pi,pi,t));
[xx0,yy0,zz0] = sph2cart(th0, phi0, radius0); % Center at (0,0,0)

% The first singluar vector
[th, phi] = meshgrid(linspace(0,2*pi,t), pi/2-theta1/180*pi);
[xx,yy,zz] = sph2cart(th, phi, radius0); % Center at (0,0,0)
[xxp,yyp,zzp] = sph2cart(0, pi/2-theta10/180*pi, radius0(1));

figure; axis equal; axis off; hold on
surface(xx0,yy0,zz0,'FaceColor','flat','FaceAlpha',0.2,'FaceLighting','gouraud','EdgeColor','none'); colormap(gray)
surface(xx,yy,zz,'FaceColor',0.3*[1 1 1],'FaceAlpha',0.2,'FaceLighting','gouraud','EdgeColor','none');

arrow3([0,0,0],[0,0,1],'q-4',3,5,1);
arrow3([0,0,0],[xxp,yyp,zzp],'z-4',3,5,1);

hold off
view([1 1 0.75])
% figurewrite('%da-rmt2',1,{0 [1 300]},figuredir);

% The second singluar vector
[th, phi] = meshgrid(linspace(0,2*pi,t), pi/2-theta2/180*pi);
[xx,yy,zz] = sph2cart(th, phi, radius0); % Center at (0,0,0)
[xxp,yyp,zzp] = sph2cart(0, pi/2-theta20/180*pi, radius0(1));

figure; axis equal; axis off; hold on
surface(xx0,yy0,zz0,'FaceColor',0.7*[1 1 1],'FaceAlpha',0.2,'FaceLighting','gouraud','EdgeColor','none'); colormap(gray)
surface(xx,yy,zz,'FaceColor','interp','FaceAlpha',0.2,'FaceLighting','gouraud','EdgeColor','none');

arrow3([0,0,0],[0,0,1],'q-4',3,5,1)
arrow3([0,0,0],[xxp,yyp,zzp],'z-4',3,5,1)

hold off
view([1 1 0.5])

% Another way to draw a sphere
% [xx, yy, zz] = meshgrid(-100:100,-100:100,-100:100);
% radius = 36;
% sphereVoxels = yy.^2 + xx.^2 + zz.^2 <= radius.^2;
% fv = isosurface(sphereVoxels, 0);
% figure;
% patch(fv,'FaceColor',0.8*[1 1 1],'FaceAlpha',0.3, 'EdgeColor','none');
% view(45,45); axis equal


% figurewrite('%da-rmt',1,{0 [1 300]},figuredir);




