% function [X, R, sigma] = ssvd(Y, svsMethod, rankEstMethod, varargin) 
% 
% Inputs:
% <Y> is the noisy matrix, could be 2D, 3D, or 4D, but assume the last
%   dimention is the repetitive measurement
% <svsMethod> determines how the singular values should be manipulated.
%   Available options are shrinkage 1 ('svs1', ref: Gavish et al 2017.
%   default), shrinkage 2 ('svs2', ref:Nadakuditi et al 2014), and truncation
%   ('none')
% <rankEstMethod> determines the method used to estimate rank and noise.
%   Available options are 'ssvd'(default), 'gavish' (ref: Gavish et al
%   2017), 'mppca' (Veraart et al 2016), 'stein' (ref: Moeller et al 2020),
%   'cordero' (Cordero-Grande et al 2019), 'sure' (ref: Ramani et al,
%   2008), 'known' (if the user provide the rank and noise information)
%
% Outputs:
% <X> is the denoised matrix
% <R> is the estimated rank
% <sigma> is the estimated noise standard deviation
%
% By Wei Zhu, zhuwei@umn.edu
%
% History:
% 2020/03/18: created this function
%
function [X, R, sigma] = ssvd(Y,svsMethod, rankEstMethod, varargin)

if ~exist('svsMethod','var') || isempty(svsMethod)
  svsMethod = 'svs1';
end
if ~exist('rankEstMethod','var') || isempty(rankEstMethod)
  rankEstMethod = 'ssvd';
end

nY = size(Y);
Y = squish(Y,length(nY)-1); % Combine the first few dimentions except the last one

% Check matrix size, if m>n, transpose
if size(Y,1) > size(Y,2)
    Y = Y.'; % conjugate transpose the matrix
    tflag = 1;
end

% Obtain matrix size
[m, n] = size(Y); % Here m<=n

% Matrix aspect ratio
beta0 = m/n;

% SVD of the input matrix
[u,sv,v] = svd(Y,'econ');

% Natural scaling
sv = diag(sv)/sqrt(n);

% ------------------- Rank and Noise Estimation ---------------------------
switch  rankEstMethod
    case 'ssvd'
        if length(varargin)<1 ||isempty(varargin{1}) % k
            k = 2;
        else
            k = varargin{1};
        end
        if length(varargin)<2 || isempty(varargin{2}) % r0
            r0 = m/2;
        else
            r0 = varargin{2};
        end
        
        % Calculate moments of the denosity distribution
        sigma0 = 1;
        qcMoment = @(x,b,k) x.^k.*sqrt((x.^2-(sigma0*(1-sqrt(b))).^2) .* ((sigma0*(1+sqrt(b))).^2-x.^2))./(pi*x*b*sigma0^2);
        moment1 = cellfun(@(k) integral(@(x) qcMoment(x,beta0,k), sigma0*(1-sqrt(beta0)), sigma0*(1+sqrt(beta0))), num2cell(k))./sigma0.^k;

        % singluar value to the kth power
        l = sv.^k; 
        csum = cumsum(l(m:-1:1,:)); % cumulative sum of "l" from small to large values, or noise --> signal for each k

        % Calculate sigma based on equation (9) in the paper
        sigma1 = (csum(m:-1:1,:)./repmat((m-(0:m-1))',1,length(k))./repmat(moment1,m,1)).^repmat((1./k),m,1); % from all component accumulation to pure noise component

        % aspect ratio here only consider the noise matrix size, which is a
        % function of variable rank r. As r increases, beta decreases and sigma
        % increases.(should we use beta0 here???? test for high rank later)
        beta = (m-(0:m-1))./(n-(0:m-1));
        % sigma2 = repmat((sv(1:m)-sv(m))./(2*sqrt(beta(:))),1,length(k)); % From singular value
        sigma2 = ((l(1:m,:)-l(m,:))./((1+sqrt(beta(:))).^k-(1-sqrt(beta(:))).^k)).^(1./k);

         %figure; plot(1:length(sigma1),sigma1,'r',1:length(sigma1),sigma2,':b','linewidth',2)

        % sigma1 > sigma2 if no signal components present in singular values;
        % sigma1 <= sigma2 if signal components present 
        rk = cellfun(@(x,y) find(x>y,1)-1, num2cell(sigma1,1), num2cell(sigma2,1),'uniformoutput',0);
        sigmak = cellfun(@(x,y) x(y+1), num2cell(sigma1,1), rk,'uniformoutput',0); % use sigma1 to estimate sigma0 is more accurate

        rk = [rk{:}];
        sigmak = [sigmak{:}];
        
%          figure; plot(rk)
%          figure; plot(sigmak)

        % Check the estimated rank and deal with extreme situations (zero matrix)
        if isempty(rk) || min(rk)>r0 % if matrix Y is close to zero matrix
            R = m;
            sigma = 0;
            % X = Y;
        else
            rk(rk>r0)=[];
            sigmak(rk>r0)=[];

            R = max(rk); % Estimated rank
            % sigma = max(sigmak(rk==R)); % Estimated noise, should close to sigma0
            sigma = max(sigmak);
            % R = rk(2); % Estimated rank
            % sigma = mean(sigmak(2)); % Estimated noise, should close to sigma0
        end
        
        
    case 'gavish'
        MPmedian = MedianMarcenkoPastur(beta0);
        sigma = median(sv) / sqrt(MPmedian);
        R = sum(sv>(1+sqrt(beta0))*sigma);
        
    case 'mppca'
        l = sv.^2; 
        csum = cumsum(l(m:-1:1,:)); % cumulative sum of "l" from small to large values, or noise --> signal for each k

        % Calculate sigmasq1 based on empirical probability distribution
        sigmasq1 = csum(m:-1:1,:)./(m-(0:m-1))'; % from all component accumulation to pure noise component
        
        % Calculate sigmasq2 based on the edges of MP distribution
        beta = (m-(0:m-1))./(n-(0:m-1));
        sigmasq2 = (l(1:m,:)-l(m,:))./((1+sqrt(beta(:))).^2-(1-sqrt(beta(:))).^2);
        
        R = find(sigmasq1>sigmasq2, 1)-1;

        if isempty(R)  % if matrix Y is close to zero matrix
            R = m;
            sigma = 0;
            % X = Y;
        else
            sigma = sqrt(sigmasq2(R+1)); % Estimated noise, should close to sigma0
        end
        
    case 'stein'
        if length(varargin)<2 ||isempty(varargin{2}) % 
            noiseMatrix = randn(5*m,5*n); % Measured noise matrix
        else
            noiseMatrix = varargin{2};
        end
        if length(varargin)<3 ||isempty(varargin{3}) % 
            nt = 10;
        else
            nt = varargin{3};
        end
        
        s1_all = zeros(nt,1);
        sigma_all = zeros(nt,1);
        for ii = 1:nt
            % Randomly sampling the noise matrix
            rowIdx = randperm(size(noiseMatrix,1), m);
            colIdx = randperm(size(noiseMatrix,2), n);
            epsilon = noiseMatrix(rowIdx,colIdx);

            stmp = svd(epsilon/sqrt(n),'econ');
            
            sigma_all(ii) = std(epsilon(:));
            s1_all(ii) = stmp(1); % Only get the largest singular value
        end
        s1 = mean(s1_all);
        
        % Obtain rank
        R = sum(sv>s1);
        % sigma = s1/(1 + sqrt(beta0));
        sigma = mean(sigma_all);
        
    case 'cordero'
        if length(varargin)<2 ||isempty(varargin{2}) % 
            sigma0 = 1; % Measured noise std
        else
            sigma0 = varargin{2};
        end
        if length(varargin)<3 ||isempty(varargin{3}) % 
            B = fix(1000/m);
        else
            B = varargin{3};
        end
        
        % Simulated noise matrix with N(0,sigma0)
        epsilon = sigma0*randn(B*m,B*n);

        stmp = svd(epsilon/sqrt(B*n),'econ');
        s1 = stmp(1);
        
        % Obtain rank
        R = sum(sv>s1);
        sigma = sigma0;
        
    case 'sure'
        sigma = varargin{2}; % Assuem known
        
        bprime = randn(m,n); % random perturbation
        epsilon = 0.00001; % scale of the perturbation bprime
        
        sure = zeros(1,m);
        for r = 1:m
            f = @(x)ssvd(x,'none','known',[r,sigma]); 
            sure(r) = MCSure(Y,f,epsilon,sigma,bprime);
        end
        % figure; plot(sure)
        
        [~,R] = min(sure);
        
    case 'known'
        if length(varargin)<2 || isempty(varargin{2})
            error('Please give the rank and noise information')
        end

        R = varargin{2}(1);
        sigma = varargin{2}(2);
        

end


% ------------------ Singular value shrinkage -------------------------
switch svsMethod
    case 'none' % Truncation
        s_opt = [sv(1:R); zeros(m-R,1)];
        
    case 'svs1'
        % Estimated matrix aspect ratio beta_hat
        if R/m <=0.05 % empirical value, what is the best?
            beta_hat = beta0;
        else
            beta_hat = (m-R)/(n-R);
        end

        beta_pos = (1+sqrt(beta_hat))*sigma; 
        beta_neg = (1-sqrt(beta_hat))*sigma;

        opt_fro_shrink = @(y)(sqrt((y.^2-beta_pos.^2).*(y.^2-beta_neg.^2))./y);
        s_opt = opt_fro_shrink(sv);
        s_opt(sv<=beta_pos) = 0;
    case 'svs2'
        r_hat = R;
        sv_noise = cat(1,diag(sv(r_hat+1:end)), zeros(n-m,m-r_hat)); % (N-r_hat) x (M-r_hat)
        w_opt = zeros(r_hat,1);
        for ii = 1:r_hat   
            sigma_i = sv(ii);
            a = (sigma_i^2*eye(n-r_hat)-sv_noise*sv_noise').^(-1);
            b = (sigma_i^2*eye(m-r_hat)-sv_noise'*sv_noise).^(-1);
            D_hat = 1/(n-r_hat)*trace(sigma_i*a)*...
                    1/(m-r_hat)*trace(sigma_i*b);

            D_hatp = 1/(n-r_hat)*trace(sigma_i*a)*...
                    1/(m-r_hat)*trace(-2*sigma_i^2 * b.^2 + b) + ...
                    1/(n-r_hat)*trace(-2*sigma_i^2 * a.^2 + a) *...
                    1/(m-r_hat)*trace(sigma_i*b);

            w_opt(ii) = -2*D_hat/D_hatp; 
        end
        s_opt = [w_opt; zeros(m-r_hat,1)];
end


%figure; plot(1:length(s),s,'ro',1:length(s),s_opt,'bs','linewidth',2)

% Signal reconstruction
X = u * diag(sqrt(n)*s_opt) * v'; 

if exist('tflag','var') && tflag == 1
    X = X.';
end
X = reshape(X, nY);

end

% --------------------------------------------------------------------

function med = MedianMarcenkoPastur(beta)
    MarPas = @(x) 1-incMarPas(x,beta,0);
    lobnd = (1 - sqrt(beta))^2;
    hibnd = (1 + sqrt(beta))^2;
    change = 1;
    while change & (hibnd - lobnd > .001),
      change = 0;
      x = linspace(lobnd,hibnd,5);
      for i=1:length(x),
          y(i) = MarPas(x(i));
      end
      if any(y < 0.5),
         lobnd = max(x(y < 0.5));
         change = 1;
      end
      if any(y > 0.5),
         hibnd = min(x(y > 0.5));
         change = 1;
      end
    end
    med = (hibnd+lobnd)./2;
end

function I = incMarPas(x0,beta,gamma)
    if beta > 1,
        error('betaBeyond');
    end
    topSpec = (1 + sqrt(beta))^2;
    botSpec = (1 - sqrt(beta))^2;
    MarPas = @(x) IfElse((topSpec-x).*(x-botSpec) >0, ...
                         sqrt((topSpec-x).*(x-botSpec))./(beta.* x)./(2 .* pi), ...
                         0);
    if gamma ~= 0,
       fun = @(x) (x.^gamma .* MarPas(x));
    else
       fun = @(x) MarPas(x);
    end
    I = quadl(fun,x0,topSpec);
    
    function y=IfElse(Q,point,counterPoint)
        y = point;
        if any(~Q),
            if length(counterPoint) == 1,
                counterPoint = ones(size(Q)).*counterPoint;
            end
            y(~Q) = counterPoint(~Q);
        end
        
    end
end

