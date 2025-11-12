
% set sampling time
T = 0.4;

s = tf('s');

amplitude = 0.1500000000;


%LAB 3 PLANT
P = -0.2588729/(s*s);


% LAB 3 (0.1s sampling time) TA plant
  % -0.001294 z - 0.001294
  % ----------------------
  %     z^2 - 2z + 1


% LAB 3 (0.4s sampling time) TA plant
  % -0.02071 z - 0.02071
  % --------------------
  %    z^2 - 2 z + 1



%Convert to Discrete
G_c2d = c2d(P, T);

% does the plant have a double integrator?
double_integrator_flag = 1;

% should the controller have an integrator?
controller_integrator_flag = 0;

%% Plant Poles and Coefficients in its Partial Fraction Decomposition

stableRealPlantPoles = [];
stableComplexPlantPoles = [];
unstablePlantPoles = [1];

if double_integrator_flag
    if unstablePlantPoles(end) ~= 1
        disp('The final unstable plant pole must be z=1!');
        stop
    elseif length(find(unstablePlantPoles == 1)) > 1
  disp('There should only be one pole at z=1 included in unstablePlantPoles!');
        stop
    end
end

stablePlantPoles = [stableRealPlantPoles stableComplexPlantPoles];
qs = [stablePlantPoles unstablePlantPoles];

%(0.4s sampling time) TA plant

% coefficents go in order of the poles
cs = [-0.02071];

if double_integrator_flag
    % coefficients include both c_n for 1/(z-1) and c_(n+1) for 1/(z-1)^2 for
    %       the pole at z=1
    c_double_integrator = -0.04142;
    cs = [cs c_double_integrator];
end     

n = length(qs);
nhat = length(stablePlantPoles);
nreal = length(stableRealPlantPoles);
ncomplex = length(stableComplexPlantPoles);

% verify that your plant is correct!
z = tf('z',T);
G = 0;
for k=1:n
    G = G + cs(k)/(z-qs(k));
end
if double_integrator_flag
    G = G + c_double_integrator/(z-1)^2;
end
G

%% Poles Chosen in the Simple Pole Approximation of W[z]

j = sqrt(-1);
realWPoles = [];

% complexWPoles = [0.4+0.1*j 0.4-0.1*j 0.5+0.1*j 0.5-0.1*j];
% % for checking the integrator in the controller:
% complexWPoles = [0.4+0.1*j 0.4-0.1*j 0.5+0.1*j 0.5-0.1*j 0.6+0.1*j 0.6-0.1*j];

%(0.4s sampling time) generate_poles(12, 0.7, 0.1)
complexWPoles = [ ...
   -0.162844740719926 - 0.112157518439654*j,  -0.162844740719926 + 0.112157518439654*j, ...
    0.219767988315681 - 0.385990883711193*j,   0.219767988315681 + 0.385990883711193*j, ...
    0.589906607603367 - 0.070650660482126*j,   0.589906607603367 + 0.070650660482126*j, ...
    0.495473959643374 + 0.412634237443595*j,   0.495473959643374 - 0.412634237443595*j, ...
    0.053602142248959 + 0.637322973169371*j,   0.053602142248959 - 0.637322973169371*j, ...
   -0.416189665998437 + 0.472808871232786*j,  -0.416189665998437 - 0.472808871232786*j ...
];



%complexWPoles = [gen_poles];

ps = [realWPoles complexWPoles];

mreal = length(realWPoles);
mcomplex = length(complexWPoles);
m = length(ps);

%% Calculation of alpha, beta, gamma, and gamma hat

alpha = zeros(m);

for i=1:m
    for k=1:n
        alpha(i,i) = alpha(i,i) + cs(k)/(ps(i)-qs(k));
    end
    if double_integrator_flag
        alpha(i,i) = alpha(i,i) + cs(n+1)/((ps(i)-1)^2);
    end
end

beta = zeros(n,m);
if double_integrator_flag
    beta = zeros(n+1,m);
end

for i=1:m
    for k=1:n
        beta(k,i) = cs(k)/(qs(k)-ps(i));
    end
    if double_integrator_flag
        beta(n,i) = beta(n,i) - cs(n+1)/((1-ps(i))^2);
        beta(n+1,i) = cs(n+1)/(1-ps(i));
    end
end

gamma = zeros(n-nhat,m);
if double_integrator_flag
    gamma = zeros(n+1-nhat,m);
end

for i=1:m
    for j=(nhat+1):n
        gamma(j-nhat,i) = cs(j)/(qs(j)-ps(i));
    end
    if double_integrator_flag
        gamma(n-nhat,i) = gamma(n-nhat,i) - cs(n+1)/((1-ps(i))^2);
        gamma(n+1-nhat,i) = cs(n+1)/(1-ps(i));
    end
end

gammaHat = zeros(n-nhat,nhat);
if double_integrator_flag
    gammaHat = zeros(n+1-nhat,nhat);
end

for k=1:nhat
    for j=(nhat+1):n
        gammaHat(j-nhat,k) = cs(j)/(qs(j)-qs(k));
    end
    if double_integrator_flag
        gammaHat(n-nhat,k) = gammaHat(n-nhat,k) - cs(n+1)/((1-qs(k))^2);
        gammaHat(n+1-nhat,k) = cs(n+1)/(1-qs(k));
    end
end

% verify on a simple example that alpha, beta, gamma, and gammahat are correct!
%alpha
%beta
%gamma
%gammaHat

%% Determination of A and b matrices for IOP equations

A = [alpha eye(m) zeros(m,nhat);
     beta [zeros(nhat,m) eye(nhat);
           zeros(size(beta,1)-nhat,m+nhat)];
     zeros(size(gamma)) gamma gammaHat];

b = [zeros(m+size(beta,1),1);
     -cs((nhat+1):end)'];

%% Determination of step response matrices

% time horizon
K = 100;

step_ry = zeros(K,m+nhat);

for k=1:K
    for i=1:m
        step_ry(k,i) = -(1-ps(i)^k)/(1-ps(i));
    end
    for j=1:nhat
        step_ry(k,m+j) = -(1-qs(j)^k)/(1-qs(j));
    end
end

step_ru = zeros(K,m);

for k=1:K
    for i=1:m
        step_ru(k,i) = (1-ps(i)^k)/(1-ps(i));
    end
end

step_ry = step_ry * amplitude;
step_ru = step_ru * amplitude;

% verify on a simple example that step_ry and step_ru are correct!
%step_ry
%step_ru

%% Determination of steady state vector

steadyState = zeros(1,m+nhat);
if controller_integrator_flag
    steadyState = zeros(3,m+nhat);
end

for i=1:m
    if ~controller_integrator_flag    
        steadyState(i) = 1/(1-ps(i));
    else
        steadyState(1,i) = 1/(1-ps(i));
        steadyState(2,i) = 1/(1-ps(i))^2;
        steadyState(3,i) = 1/(1-ps(i))^3;
    end
end

for k=1:nhat
    if ~controller_integrator_flag
        steadyState(m+k) = 1/(1-qs(k));
    else
        steadyState(1,m+k) = 1/(1-qs(k));
        steadyState(2,m+k) = 1/(1-qs(k))^2;
        steadyState(3,m+k) = 1/(1-qs(k))^3;
    end
end

% verify on a simple example that steadyState is correct!
%steadyState

%% Defining the variables for the optimization

wreal = sdpvar(mreal,1,'full');
wcomplex = sdpvar(mcomplex/2,1,'full','complex');
w = wreal;
for i=1:(mcomplex/2)
    w = [w;
         wcomplex(i);
         conj(wcomplex(i))];
end

xreal = sdpvar(mreal,1,'full');
xcomplex = sdpvar(mcomplex/2,1,'full','complex');
x = xreal;
for i=1:(mcomplex/2)
    x = [x;
         xcomplex(i);
         conj(xcomplex(i))];
end

xhatreal = sdpvar(nreal,1,'full');
xhatcomplex = sdpvar(ncomplex/2,1,'full','complex');
xhat = xhatreal;
for i=1:(ncomplex/2)
    xhat = [xhat;
            xhatcomplex(i);
            conj(xhatcomplex(i))];
end


%% Defining the objective function and constraints for the optimization

Objective = 0;

% IOP constraint
Constraints = [A*[w;x;xhat] == b];

% input saturation constraint
Constraints = [Constraints,
               max(step_ru*w) <= 0.7,
               min(step_ru*w) >= -0.7];

% steady state constraint
if ~controller_integrator_flag
    Constraints = [Constraints,
                   steadyState*[x;xhat]+1==0];
else
    Constraints = [Constraints,
                   steadyState*[x;xhat]+[1;0;0]==[0;0;0]];
end

% overshoot constraint
Constraints = [Constraints,
               max(step_ry*[x;xhat]) <= amplitude*1.45*(-steadyState(1,:)*[x;xhat])];

% settling time constraint
jhat = 7/T;
Constraints = [Constraints,
               max(step_ry(jhat:end,:)*[x;xhat]) <= amplitude*1.02*(-steadyState(1,:)*[x;xhat]),
               min(step_ry(jhat:end,:)*[x;xhat]) >= amplitude*0.98*(-steadyState(1,:)*[x;xhat])];

%% Solving the optimization problem

% set some options for YALMIP and solver
options = sdpsettings('verbose',1,'solver','mosek');

% solve the problem
sol = optimize(Constraints,Objective,options);

% obtain the solution
wsol = value(w);
xsol = value(x);
xhatsol = value(xhat);

%% Plotting the solution

figure(1)
plot(T*(1:K),step_ry*[xsol;xhatsol]);
xlabel('Time [s]');
ylabel('y[k]');

figure(2)
plot(T*(1:K),step_ru*wsol);
xlabel('Time [s]');
ylabel('u[k]');

% use log scale for heat map?
log_scale_flag = 1;

% heat map
figure(3)
t = linspace(0,2*pi);
plot(cos(t),sin(t),'k--');
hold on;
if log_scale_flag
    scatter(real(ps),imag(ps),50,log(abs(wsol)),'filled');
    scatter(real(qs(1:nhat)),imag(qs(1:nhat)),50,log(abs(xhatsol)),'filled');
else
    scatter(real(ps),imag(ps),50,abs(wsol),'filled');
    scatter(real(qs(1:nhat)),imag(qs(1:nhat)),50,abs(xhatsol),'filled');
end
hold off;
colormap(jet);
colorbar;

%% Recover the transfer functions

z = tf('z',T);

% calculate W
W = 0;
for i=1:m
    W = W + wsol(i)/(z-ps(i));
end

% calculate X
X = 1;
for i=1:m
    X = X + xsol(i)/(z-ps(i));
end
for k=1:nhat
    X = X + xhatsol(k)/(z-qs(k));
end

% remove the imaginary coefficients in W
[num,den] = tfdata(W);
num{1} = real(num{1});
den{1} = real(den{1});
W = tf(num,den,T);

% remove the imaginary coefficients in X
[num,den] = tfdata(X);
num{1} = real(num{1});
den{1} = real(den{1});
X = tf(num,den,T);

% find the poles and zeros of W and X (if desired)
%zpk(W)
%zero(W)
%pole(W)
%zpk(X)
%zero(X)
%pole(X)

%% Calculate D[z]
format long
D = W / X;

%Show simplified D[z]
zpk(D)

[num_d, den_d] = tfdata(D,'v');  




%% Verify design in DT

% % compute D by hand
% j = sqrt(-1);
% D = (0.15246*(z-0.8423)*(z-0.8))/((z+0.4903)*(z-0.7796)*(z-0.3107));
% 
% % compute T_ry and T_ru by hand  (using Nf, Dg, etc)
% T_ry = (0.15246*(z-0.8423)*(z-0.8)*5*(z-0.7672)*(z-0.3128))/...
%        (0.15246*(z-0.8423)*(z-0.8)*5*(z-0.7672)*(z-0.3128) + ...
%         (z+0.4903)*(z-0.7796)*(z-0.3107)*(z-1)^2*(z-0.8));
% T_ru = (0.15246*(z-0.8423)*(z-0.8)*(z-1)^2*(z-0.8))/...
%        (0.15246*(z-0.8423)*(z-0.8)*5*(z-0.7672)*(z-0.3128) + ...
%         (z+0.4903)*(z-0.7796)*(z-0.3107)*(z-1)^2*(z-0.8));
% 
% figure(1)
% hold on;
% step(T_ry,'g');
% hold off;
% 
% figure(2)
% hold on;
% step(T_ru,'g');
% hold off;


%% Metrics for the step response (using solved w/x/xhat)

% model-predicted trajectories (your step_ry/step_ru are already scaled by 'amplitude')
y_traj = real(step_ry*[xsol; xhatsol]);   % Kx1
u_traj = real(step_ru*wsol);              % Kx1

% step endpoints
y_final_desired = amplitude;              % desired final value for the step
y_ss            = y_traj(end);            % achieved steady-state (last point in horizon)

% steady-state error
ess = y_ss - y_final_desired;

% overshoot
[y_max, k_peak] = max(y_traj);
OS = (y_max - y_final_desired) / max(1e-12, amplitude);   % safe divide

% 2% settling time (relative to desired final value)
tol = 0.02 * amplitude;
Ksettle = K;  % default: assume settled by end
for kk = K-1:-1:1
    if abs(y_traj(kk) - y_final_desired) > tol
        Ksettle = kk + 1;   % first index from which we're within band
        break;
    end
end
Ts = T * (Ksettle - 1);     % time at sample Ksettle

% input usage
maxu = max(abs(u_traj));

% Pretty print
fprintf('\n---- Metrics (0 → %.4g) ----\n', amplitude);
fprintf('Desired final value       : %.6g\n', y_final_desired);
fprintf('Achieved steady state y_ss: %.6g\n', y_ss);
fprintf('Steady-state error (ess)  : %.6g\n', ess);
fprintf('Peak value y_max (k=%d)   : %.6g\n', k_peak, y_max);
fprintf('Overshoot                 : %.3f %%\n', 100*max(0,OS));
fprintf('2%% settling index         : %d of %d\n', Ksettle, K);
fprintf('2%% settling time Ts       : %.6g s\n', Ts);
fprintf('Max |u[k]|                : %.6g\n', maxu);


