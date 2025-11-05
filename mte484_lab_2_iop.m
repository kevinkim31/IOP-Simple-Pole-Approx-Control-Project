%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This template is for the IOP control design with SPA, and is incomplete.
% You need to complete it by replacing every * with the correct code.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set time step
T = 0.02;
j = sqrt(-1);

amplitude = 1.4000;
offset = -0.7;

s = tf('s');
P = 2.51124/(s*(0.02821*s+1)); % own values
%P = 2.542/(s*(0.022*s+1)); % TA values


%--- Pick Ts from rule of thumb: ws > (5..10)w_bw ---
% poles = pole(P);                  % continuous-time poles (rad/s)
% w_bw  = max(abs(real(poles)));    % bandwidth proxy from fastest decay
% fac   = 10;                       % choose 5..10
% ws    = fac * w_bw;               % sampling ang. freq (rad/s)
% SF    = 3.5;
% T     = 2*pi/ws/SF;                  % sampling period (s)



G_c2d = c2d(P, T);

syms z

% with 0.02 ms sampling time
  %  0.01425 z + 0.01126
  % ----------------------
  % z^2 - 1.492 z + 0.4922

% with 0.015 ms sampling time
  % 0.008453 z + 0.007082
  % ----------------------
  % z^2 - 1.588 z + 0.5876

% with 0.01ms sampling time
  % 0.003969 z + 0.003527
  % ----------------------
  % z^2 - 1.702 z + 0.7015

% with 0.01ms sampling time (TA plant)
  % 0.004993 z + 0.004292
  % ----------------------
  % z^2 - 1.635 z + 0.6347

% with SF = 1.1 (0.0161ms sampling time)
  % 0.009638 z + 0.007971
  % ----------------------
  % z^2 - 1.565 z + 0.5648

% with 0.006ms sampling time
  % 0.001495 z + 0.001392
  % ----------------------
  % z^2 - 1.808 z + 0.8084

% with 0.008ms sampling time
  % 0.002597 z + 0.002363
  % ----------------------
  % z^2 - 1.753 z + 0.7531

% with 0.006ms sampling time (TA plant)
  % 0.001903 z + 0.001738
  % ----------------------
  % z^2 - 1.761 z + 0.7613

%% Plant Poles and Coefficients in its Partial Fraction Decomposition

% test plant (T = 0.02):
stableRealPlantPoles = [0.999606 0.492394];
stableComplexPlantPoles = [];
unstablePlantPoles = [];
cs = [0.0502835 -0.0360335]; %coefficients

% test plant (T = 0.015):
% stableRealPlantPoles = [0.5870314033];
% stableComplexPlantPoles = [];
% unstablePlantPoles = [1.0009685967];
% cs = [-0.0290966278 0.0375496278]; %coefficients

%test plant (T = 0.01):
% stableRealPlantPoles = [0.700331];
% stableComplexPlantPoles = [];
% unstablePlantPoles = [1.00167];
% cs = [-0.0209288 0.0248978 ]; %coefficients

% %test plant (T = 0.01) (TA plant):
% stableRealPlantPoles = [0.63418];
% stableComplexPlantPoles = [];
% unstablePlantPoles = [1.00082];
% cs = [-0.0203427 0.0253357]; %coefficients

% %test plant (SF = 1.1; T = 0.0161)
% stableRealPlantPoles = [0.564541];
% stableComplexPlantPoles = [];
% unstablePlantPoles = [1.00046];
% cs = [-0.0307673 0.0404053]; %coefficients

% %test plant (T = 0.006):
% stableRealPlantPoles = [0.997894 0.810106];
% stableComplexPlantPoles = [];
% unstablePlantPoles = [];
% cs = [0.015357 -0.013862]; %coefficients

% %test plant (T = 0.006) (TA plant):
% stableRealPlantPoles = [0.998738 0.762262];
% stableComplexPlantPoles = [];
% unstablePlantPoles = [];
% cs = [0.0153867 -0.0134837]; %coefficients

%test plant (T = 0.008):
% stableRealPlantPoles = [0.999594 0.753406];
% stableComplexPlantPoles = [];
% unstablePlantPoles = [];
% cs = [0.0201428 -0.0175458]; %coefficients

stablePlantPoles = [stableRealPlantPoles stableComplexPlantPoles];
qs = [stablePlantPoles unstablePlantPoles];
qs

cs

n = length(qs);
nhat = length(stablePlantPoles);
nreal = length(stableRealPlantPoles);
ncomplex = length(stableComplexPlantPoles);


%% Poles Chosen in the Simple Pole Approximation of W[z]

realWPoles = [];

%default (6 poles)
%complexWPoles = [-0.229587+0.0979663*j -0.229587-0.0979663*j 0.104614+0.337152*j 0.104614-0.337152*j 0.427919+0.0617113*j 0.427919-0.0617113*j]; % only 6 closest?

%T = 0.02 (8 poles) 0.8 range
%complexWPoles =[-0.367905997907408+0.156987823425114*j -0.367905997907408-0.156987823425114*j 0.167640262224593+0.540274691690505*j 0.167640262224593-0.540274691690505*j 0.685726405854943+0.098890324669614*j 0.685726405854943-0.098890324669614*j 0.553548232962454+0.577567618365283*j 0.553548232962454-0.577567618365283*j ]

%T = 0.02 (6 poles) 0.3 range
complexWPoles =[ -0.159307970196240-0.067977721585487*j -0.159307970196240+0.067977721585487*j 0.072590362891791-0.233945804012891*j 0.072590362891791+0.233945804012891*j 0.296928243758090-0.042820766676188*j 0.296928243758090+0.042820766676188*j]


%T = 0.008 (6 poles) 0.9 range
%complexWPoles = [-0.213734+0.0912017*j -0.213734-0.0912017*j 0.0973902+0.313871*j 0.0973902-0.313871*j 0.398371+0.0574501*j 0.398371-0.0574501*j]


%from gen poles function
%complexWPoles = [gen_poles];


ps = [realWPoles complexWPoles];

ps

mreal = length(realWPoles);
mcomplex = length(complexWPoles);
m = length(ps);

%% Calculation of alpha, beta, gamma, and gamma hat

alpha = zeros(m);

for i=1:m
    for k=1:n
        alpha(i,i) = alpha(i,i) + cs(k)/(ps(i)-qs(k));
    end
end

beta = zeros(n,m);

for i=1:m
    for k=1:n
        beta(k,i) = cs(k)/(qs(k)-ps(i));
    end
end

gamma = zeros(n-nhat,m);

for i=1:m
    for j=(nhat+1):n
        gamma(j-nhat,i) = cs(j)/(qs(j)-ps(i));
    end
end

gammaHat = zeros(n-nhat,nhat);

for k=1:nhat
    for j=(nhat+1):n
        gammaHat(j-nhat,k) = cs(j)/(qs(j)-qs(k));
    end
end

% TODO: verify on a simple example that alpha, beta, gamma, and gammahat are correct!
alpha
beta
gamma
gammaHat

%% Determination of A and b matrices for IOP equations

A = [alpha eye(m) zeros(m,nhat);
     beta [zeros(nhat,m) eye(nhat);
           zeros(size(beta,1)-nhat,m+nhat)];
     zeros(size(gamma)) gamma gammaHat];

b = [zeros(m+size(beta,1),1);
     -cs((nhat+1):end) ];

A
b

%% Determination of step response matrices

% time horizon
K = 200;
% K=5; % initial K

amplitude = 1.4;
offset = -0.7;

step_ry = zeros(K,m+nhat);

for k=1:K
    for i=1:m
        step_ry(k,i) = -(1-ps(i)^k)/(1-ps(i)); % row k, i=1 to m
    end
    for j=1:nhat
        step_ry(k,m+j) = -(1-qs(j)^k)/(1-qs(j));
    end
end

% step_ry=amplitude*step_ry;

step_ru = zeros(K,m);

for k=1:K
    for i=1:m
        step_ru(k,i) = (1-ps(i)^k)/(1-ps(i));
    end
end


step_ry = step_ry * amplitude;
step_ru = step_ru * amplitude;

% step_ru=amplitude*step_ru;

% verify on a simple example that step_ry and step_ru are correct!
% step_ry
% step_ru

%% Determination of steady state vector

steadyState = zeros(1,m+nhat);

for i=1:m
    steadyState(i) = 1/(1-ps(i));
end

for k=1:nhat
    steadyState(m+k) = 1/(1-qs(k));
end

% verify on a simple example that steadyState is correct!
% steadyState = amplitude*steadyState;

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

%Objective = 0;
Objective = 0;

fprintf('size(A)=%dx%d, size(w)=%dx%d, size(x)=%dx%d, size(xhat)=%dx%d, size(b)=%dx%d\n', ...
    size(A), size(w), size(x), size(xhat), size(b));

% IOP constraint
Constraints = [A*[w;x;xhat] == b];


% input saturation constraint
Constraints = [Constraints,
               max(step_ru*w) <= 6,
               min(step_ru*w) >= -6];

% steady state constraint
Constraints = [Constraints,
               1+steadyState*[x;xhat] == 0];

% overshoot constraint
Constraints = [Constraints,
               max(step_ry*[x;xhat]) <= (amplitude+0.05)*(-steadyState*[x;xhat])];

jhat = 0.2/T;
% settling time constraint
Constraints = [Constraints,
               max(step_ry(jhat:end, :)*[x;xhat]) <= amplitude*1.02*(-steadyState*[x;xhat]),
               min(step_ry(jhat:end, :)*[x;xhat]) >= amplitude*0.98*(-steadyState*[x;xhat])];
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
title('y[k] output')
ylabel('y[k]');

figure(2)
plot(T*(1:K),step_ru*wsol);
xlabel('Time [s]');
title('u[k] output')
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

W
X

% --- Using your existing W and X transfer functions ---
D = W / X;

% hand num and den of D[z] over to simulink discrete TF block
format long

% Show the simplified discrete transfer function
zpk(D)

[zD0, pD0, ~] = zpkdata(D, 'v');
fprintf('\nBefore simplification:\n');
fprintf('  Zeros: %d\n', numel(zD0));
fprintf('  Poles : %d\n', numel(pD0));


%% Snap nearly identical poles/zeros (e.g. z=1 vs z=0.998)
[zD_snap, pD_snap, kD_snap] = zpkdata(D, 'v');
snap_tol = 0.002;  % adjust if needed
pD_snap(abs(pD_snap - 1) < snap_tol) = 1;
zD_snap(abs(zD_snap - 1) < snap_tol) = 1;
D_snap = zpk(zD_snap, pD_snap, kD_snap, D.Ts);

%% Simplify D[z] by canceling pole-zero pairs and display result
tol = 0.001;                  % numerical tolerance
D_simplified = minreal(D_snap, tol);  % remove near-canceling poles/zeros
disp('Final simplified D[z]:');
zpk(D_simplified)

% --- Count poles and zeros of final D[z] ---
[zD, pD, kD] = zpkdata(D_simplified, 'v');
fprintf('\nAfter simplification (tol = %.3g):\n', tol);
fprintf('\nNumber of zeros in D[z]: %d\n', numel(zD));
fprintf('Number of poles in D[z]: %d\n', numel(pD));



[num_d, den_d] = tfdata(D_simplified,'v');  
% num_d = num_d/den_d(1);
% den_d = den_d/den_d(1);

% find the poles and zeros of W and X
zpk(W);
zero(W);
pole(W);
zpk(X);
zero(X);
pole(X);

%% Verify design in discrete time
% 
% % compute D by hand
% j = sqrt(-1);
% D = 0.43333*((z-0.4)*(z-0.5))/((z-1)*(z-0.3393));
% zpk(D)
% 
% denom=1+G*D;
% T_ry_num=G*D;
% zpk(T_ry_num)
% zpk(denom)
% 
% % % compute T_ry and T_ru by hand
% T_ry = 0.88399*(z-0.4)/((z+0.02002)*(z^2-0.8*z+0.32));
% T_ru = 0.88399*(z-0.4)*(z-0.5)/((z+0.02002)*(z^2-0.8*z+0.32));
% 
% 
% t = T*(1:K);
% 
% 
% figure(3)
% [y1, t1] = step(T_ry, t);
% plot(t1, y1), grid on
% xlabel('Time [s]'); ylabel('T_{ry}[k]'); title('Step: T_{ry}')
% 
% figure(4)
% [y2, t2] = step(T_ru, t);
% plot(t2, y2), grid on
% xlabel('Time [s]'); ylabel('T_{ru}[k]'); title('Step: T_{ru}')


% ---------- Values & metrics: 0 → amplitude ----------
wsol  = value(w);
xsol  = value(x);
xhsol = value(xhat);

% model-predicted trajectories (no offset)
y_traj = real(step_ry*[xsol; xhsol]);   % Kx1, step_ry already scaled by 'amplitude'
u_traj = real(step_ru*wsol);            % Kx1, step_ru already scaled by 'amplitude'

% step endpoints
y0                = 0;
y_final_desired   = amplitude;          % e.g., 1.4

% achieved steady state (last sample in the horizon)
y_ss = y_traj(end);

% steady-state error (signed, relative to desired)
ess = y_ss - y_final_desired;

% overshoot (relative to commanded amplitude)
[y_max, k_peak] = max(y_traj);
OS = (y_max - y_final_desired) / max(1e-12, amplitude);   % clip to 0 when reporting in %

% 2% settling time (band around desired final)
tol = 0.02 * amplitude;
Ksettle = K;  % default: last sample is settled
for kk = K-1:-1:1
    if abs(y_traj(kk) - y_final_desired) > tol
        Ksettle = kk + 1;               % first index from which we're within band
        break;
    end
end
Ts = T * (Ksettle - 1);                 % time at sample Ksettle


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
fprintf('Max |u[k]|                 : %.6g\n', maxu);