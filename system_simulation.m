%% System Simulation Assignment
b = .2; % Axle width (distance from left wheel to right wheel)
vl = [0 .1 .2 .3 .4 .5 .5 .5 .5 .5 .6 .7 .8 .9 .9 .9 .9 .9 .7 .5 .3 .1 0]; % Measured velocities of left wheel
vr = [0 .1 .2 .3 .4 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .3 .1 0]; % Measured velocities of right wheel
t = 0:1:22; % Time vector
% Summary of desired path:
%  Accelerate gradually forward for 5 seconds. 
%    Hold a steady speed forward for 4 seconds.
%    Turn left for 9 seconds (gradually, then at a steady turn, then
%      gradually coming back to forward motion)
%    Decelerate to a stop. 

% Visualization of wheel velocities: 
tiledlayout('flow');
nexttile;
plot(t, vl);
title('Left wheel velocity vs. time')
nexttile;
plot(t, vr);
title('Right wheel velocity vs. time')

% Create timeseries object for each wheel's velocities
tsl = timeseries(vl, t); 
tsr = timeseries(vr, t);

L = 23; % Number of samples in vl, vr
mu = 0; % Mean for noise in x and y measurements
sigma_xy = .05; % Covariance (or standard deviation?)
sigma_phi = .05; % Covariance for angle noise
x_noise = sigma_xy * randn(L,1); % add mu for non-zero mean
y_noise = sigma_xy * randn(L,1);
phi_noise = sigma_phi * randn(L,1);

tsx = timeseries(x_noise, t);
tsy = timeseries(y_noise, t);
tsphi = timeseries(phi_noise, t);

% Left disturbances
dist_l = zeros(23, 1);
dist_l(5:6) = .2;

% Right disturbances
dist_r = zeros(23,1);
dist_r(17:20) = .15;

tsdl = timeseries(dist_l, t);
tsdr = timeseries(dist_r, t);

% Note from Dr. Vincent: 
% Hi David, the simulation looks good. In addition to having a direct 
%   measurement of x and y, I'd like you to add measurements of the distance 
%   to beacons at known locations. Pick 3 fixed locations (e.g. p1, p2 and p3) 
%   and report norm([x1;x2]-p_i) for i=1 to 3. When you do the estimation, 
%   have state estimation using just the beacons as one case.

%% State Estimation with Extended Kalman Filter (EKF)
b = .2;
% vr = 1;
% vl = 1;
% Initial Conditions
Q = diag([.1 .1 .1]); % Covariance matrix for plant noise (?)
R = diag([.1 .1 .1]); % Covariance matrix for measurement noise (?)
P0 = diag([0 0 0]); % Covariance of initial state (0 since known exactly)
x0 = [0;0;0]; % Initial coordinates
Pk_plus = P0; 
xk_plus = x0;
syms x y phi
% f = [x + cos(phi)*vr/2 + cos(phi)*vl/2;
%     y + sin(phi)*vr/2 + sin(phi)*vl/2;
%     phi + vr/b - vl/b];

%% Moved into separate function file (this is outdated)
for t = 0:1:22
    %    Time Update section:

    % Plug in and calculate jacobian:
    f = [x + cos(phi)*vr(t+1)/2 + cos(phi)*vl(t+1)/2;
        y + sin(phi)*vr(t+1)/2 + sin(phi)*vl(t+1)/2;
        phi + vr(t+1)/b - vl(t+1)/b];
    bigPhi = jacobian(f,[x y phi]);
    % substitute states
    bigPhi = subs(bigPhi, [x; y; phi], xk);

    % Update state
    xkplus1_minus = subs(f, [x;y;phi], xk);

    % Update Covariance
    Pkplus1_minus = bigPhi * Pk * bigPhi' + Q;

    Pk = Pkplus1_minus;
    xk = xkplus1_minus;
    
    
    %     Measurement Update section:
    g = [x;
        y;
        phi];
    Ck = jacobian(g, [x y phi]);
    % Gives: Ck = eye(3) (for now)
    
    % Kalman gain:
    Kk = Pk * Ck.' * inv(Ck * Pk * Ck.' + R);
    xkplus = xkminus + Kk*(yk - subs(g, xk));
    Pkplus = (eye(3,3) - Kk*Ck)*Pk;

end
