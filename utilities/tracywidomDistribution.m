% Theory : Compute and plot the Tracy-Widom distribution
% 
% Reference: A. Edelman, et al, "Random matrix theory, numerical
% computation and applications". Proc Sym Ap, vol 72, pp. 53-82, 2014
%
% By Wei Zhu, zhuwei@umn.edu

function [f1, f2, f4] = tracywidomDistribution(xrange,tau)

% Parameters
if ~exist('xrange','var') || isempty(xrange)
  xrange = -8:0.005:5;
  tau = 1;
end 
t0=xrange(end); % right endpoint

% Theory : The differential equation solver
deq=@(t,y) [y(2); t*y(1)+ 2*y(1)^3; y(4); y(1)^2];
opts=odeset ('reltol' ,1e-12, 'abstol' ,1e-15);
y0=[airy(t0); airy(1,t0); 0; airy(t0)^2]; % boundary conditions
[t ,y]=ode45(deq, flip(xrange), y0, opts) ; % solve

%q = interp1(t,y(:,1),t*tau);
q = y(:,1);

dI = -[0; cumsum((q(1:end-1).^2+q(2:end).^2)/2.*diff(t))]; 
I = -[0; cumsum((dI(1:end-1)+dI(2:end))/2.*diff(t))];
J = -[0; cumsum(q(1:end-1)+q(2:end))/2.*diff(t)];

F2 = exp(-I/tau^2);
F1 = sqrt(F2.*exp(-J/tau));
F4 = sqrt(F2).*(exp(J/2/tau)+exp(-J/2/tau))/2;
t4 = t*tau/2^(2/3);

f2 = gradient(F2,t*tau);
f1 = gradient(F1,t*tau);
f4 = gradient(F4,t4*tau);

% The following results are very close to F2 and f2 obtained above
% F2=exp(-y(:,3)); % the distribution
% f2=gradient(F2, t) ; % the density

% Plot
% figure; plot(t, f2,t,f1,t4,f4, 'LineWidth', 2)
% axis([-5 2 0 .5])


