%% Example EKF Implementation assignment
%% (mostly) from System Simulation Assignment:
b = .2; % Axle width (distance from left wheel to right wheel)
piv1 = .05 * pi / 3;
piv2 = -1 * (piv1 + piv1 * .927/1.57);
vl = [1 1 1 1 1 -piv1 -piv1 -piv1 1 1 1 1 1 -piv2 -piv2 -piv2 1 1 1 1 1 0 0]; % Measured velocities of left wheel
vr = [1 1 1 1 1  piv1  piv1  piv1 1 1 1 1 1  piv2  piv2  piv2 1 1 1 1 1 0 0]; % Measured velocities of right wheel
t = 0:1:22; % Time vector
% Summary of desired path:
% Start at origin
%  Arrive at beacon one (5,0) at k = 5
%  Pivot towards beacon two (5,5) for 3 seconds
%  Arrive at beacon two at k = 13
%  Pivot towards beacon three (8,1) for 3 seconds 
%  Arrive at beacon three at k = 21

% Visualization of wheel velocities:
% tiledlayout('flow');
% nexttile;
% plot(t, vl);
% title('Left wheel velocity vs. time')
% nexttile;
% plot(t, vr);
% title('Right wheel velocity vs. time')

% Create timeseries object for each wheel's velocities
tsl = timeseries(vl, t);
tsr = timeseries(vr, t);

L = 23; % Number of samples in vl, vr

mu = 0; % Mean for process noise
sigma_proc_xy = .05; % Covariance for process noise
sigma_proc_phi = .02; % Covariance for process noise
x_proc_noise = sigma_proc_xy * randn(L,1); % If non-zero mean ever desired, add mu
y_proc_noise = sigma_proc_xy * randn(L,1);
phi_proc_noise = sigma_proc_phi * randn(L,1);

mu2 = 0; % Mean for noise in x and y measurements
sigma_meas_xy = .15; % Covariance for measurement noise
sigma_meas_phi = .05; % Covariance for measurement noise
x_meas_noise = sigma_meas_xy * randn(L,1); % If non-zero mean ever desired, add mu2
y_meas_noise = sigma_meas_xy * randn(L,1);
phi_meas_noise = sigma_meas_phi * randn(L,1);

% Create timeseries objects for process noise (simulink)
tsx = timeseries(x_noise, t);
tsy = timeseries(y_noise, t);
tsphi = timeseries(phi_noise, t);

% Create timeseries object for measurement (simulink)
tsxm = timeseries(x_meas_noise, t);
tsym = timeseries(y_meas_noise, t);
tsphim = timeseries(phi_meas_noise, t);

% Disturbance: 
dist_r = zeros(23,1);
dist_r(19:21) = -.1;

tsdr = timeseries(dist_r, t);

%% State Estimation with Extended Kalman Filter (EKF)
b = .2;
% vr = 1;
% vl = 1;
% Initial Conditions

Q = diag([sigma_proc_xy sigma_proc_xy sigma_proc_phi]); % Covariance matrix for plant/process noise
R = diag([sigma_meas_xy sigma_meas_xy sigma_meas_phi]); % Covariance matrix for measurement noise
P0 = diag([0 0 0]); % Covariance of initial state (0 since known exactly)
x0 = [0;0;0]; % Initial coordinates
Pk_plus = P0;
xk_plus = x0;
syms x y phi
% f = [x + cos(phi)*vr/2 + cos(phi)*vl/2;
%     y + sin(phi)*vr/2 + sin(phi)*vl/2;
%     phi + vr/b - vl/b];

%% 
xDot=0;
yDot=0;
phiDot=0;

xpos_ = zeros(1,23); % For noiseless sim
xpos = zeros(1,23);

ypos_ = zeros(1,23); % Noiseless sim
ypos = zeros(1,23);

ykn = [zeros(1,23); zeros(1,23); zeros(1,23)]

phiPos = zeros(1,23);
phiPos_ = zeros(1,23);

estimate = zeros(3,23);

x_tris = zeros(1,23);
y_tris = zeros(1,23);

xlast_ = 0; % for noiseless sim section
ylast_ = 0;
phiLast_ = 0;

xlast = 0;
ylast = 0;
phiLast = 0;

p1 = [5; 0] % beacon 1 location
p2 = [5; 5] % beacon 2 location
p3 = [8; 1] % beacon 3 location

% Error measurement matrices
simError = zeros(23,3);
estError = zeros(23,3);

distance = zeros(23,1);
estDistance = distance;

for i = 0:1:22
    % Ideal Simulation section (No disturbance, no noise)
    xDot_ = .5 * cos(phiLast_) * (vr(i+1) + vl(i+1));
    xlast_ = xDot_ + xlast_;
    xpos_(i+1) = xlast_; % Record in array
    
    yDot_ = .5 * sin(phiLast_)* (vr(i+1) + vl(i+1));
    ylast_ = yDot_ + ylast_;
    ypos_(i+1) = ylast_; % Record in arry
    
    phiDot_ = vr(i+1)/b - vl(i+1)/b;
    phiLast_ = phiDot_ + phiLast_;
    phiPos_(i+1) = phiLast_; % Record in array
    
    
    % Simulation section:
    xDot = .5 * cos(phiLast)* (vr(i+1) + dist_r(i+1) + vl(i+1)) + x_proc_noise(i+1);
    xlast = xDot + xlast;
    xpos(i+1) = xlast;
    
    yDot = .5 * sin(phiLast)* (vr(i+1) + dist_r(i+1) + vl(i+1)) + y_proc_noise(i+1);
    ylast = yDot + ylast;
    ypos(i+1) = ylast;
    
    phiDot = 1/b * (vr(i+1) + dist_r(i+1) - vl(i+1)) + phi_proc_noise(i+1);
    phiLast = phiDot + phiLast;
    phiPos(i+1) = phiLast;
    
    % Add measurement noise:
    yk = [xpos(i+1); ypos(i+1); phiPos(i+1)] + [x_meas_noise(i+1); y_meas_noise(i+1); phi_meas_noise(i+1);];

    
    %    Time Update section:

    % Plug in and calculate jacobian:
    % f = [x + cos(phi)*vr(i+1)/2 + cos(phi)*vl(i+1)/2;
    %      y + sin(phi)*vr(i+1)/2 + sin(phi)*vl(i+1)/2;
    %      phi + vr(i+1)/b - vl(i+1)/b];
    % bigPhi = jacobian(f, [x y phi]);
    
    % Substitute states:
    % bigPhi = subs(bigPhi, [x; y; phi], xk_plus);
    bigPhi = [1 0 (-0.5000)*vl(i+1)*sin(xk_plus(3))-0.5000*vr(i+1)*sin(xk_plus(3));
        0 1 0.5000*vl(i+1)*cos(xk_plus(3)) + 0.5000*vr(i+1)*cos(xk_plus(3));
        0 0 1];

    % Update state
    %xkplus1_minus = subs(f, [x;y;phi], xk_plus);
    xk_plus1_minus = [xk_plus(1) + cos(xk_plus(3))*vr(i+1)/2 + cos(xk_plus(3))*vl(i+1)/2;
        xk_plus(2) + sin(xk_plus(3))*vr(i+1)/2 + sin(xk_plus(3))*vl(i+1)/2;
        xk_plus(3) + vr(i+1)/b - vl(i+1)/b];

    % Update Covariance
    Pkplus1_minus = bigPhi * Pk_plus * bigPhi' + Q;

    Pk = Pkplus1_minus;
    xk_minus = xk_plus1_minus;
   
    %     Measurement Update section:
    g = [x;
        y;
        phi];
    %Ck = jacobian(g, [x y phi]);
    Ck = eye(3,3);
    % Gives: Ck = eye(3) (for now)
   
    % Kalman gain:
    Kk = Pk * Ck.' * inv(Ck * Pk * Ck.' + R);
    %xk_plus = xk_minus + Kk*(yk - subs(g,[x;y;phi], xk_minus));
    xk_plus = xk_minus + Kk*(yk - xk_minus);
    Pk_plus = (eye(3,3) - Kk*Ck)*Pk;
    estimate(:,i+1) = xk_plus;
    
    %% Beacon section
    % Distance with noisy simulated measurements
    d_1 = norm([xpos(i+1); ypos(i+1)]-p1); % Distance to beacon 1
    d_2 = norm([xpos(i+1); ypos(i+1)]-p2); % Distance to beacon 2
    d_3 = norm([xpos(i+1); ypos(i+1)]-p3); % Distance to beacon 2
    simError(i+1, :) = [d_1, d_2, d_3];
    
    % Use beacon distance measurements to triangulate robot position: 
    [xout,yout] = circcirc(p1(1), p1(2), d_1, p2(1), p2(2), d_2);
    [xout2,yout2] = circcirc(p1(1), p1(2), d_1, p3(1), p3(2), d_3);
    if norm([xout(1); yout(1)]-[xout2(1); yout2(1)]) < .01
        x_tri = xout(1);
        y_tri = yout(1);
    end
    if norm([xout(2); yout(2)]-[xout2(1); yout2(1)]) < .01
        x_tri = xout(2);
        y_tri = yout(2);
    end
    if norm([xout(1); yout(1)]-[xout2(2); yout2(2)]) < .01
        x_tri = xout(1);
        y_tri = yout(1);
    end
    if norm([xout(2); yout(2)]-[xout2(2); yout2(2)]) < .01
        x_tri = xout(2);
        y_tri = yout(2);
    end
     
    x_tris(i+1) = x_tri;
    y_tris(i+1) = y_tri;
     
%    Distance with EKF estimations
%    estd_1 = norm([xk_plus(1);xk_plus(2)]-p1); % Distance to beacon 1
%    estd_2 = norm([xk_plus(1);xk_plus(2)]-p2); % Distance to beacon 2
%    estd_3 = norm([xk_plus(1);xk_plus(2)]-p3); % Distance to beacon 2
%    estError(i+1, :) = [estd_1 estd_2 estd_3];
    
    % Distance from ideal path:
    d = norm([xpos(i+1); ypos(i+1)]-[xpos_(i+1); ypos_(i+1)]); % Distance to ideal path
    estd = norm([xk_plus(1); xk_plus(2)]-[xpos(i+1); ypos(i+1)]); % Distance to noisy path
    distance(i+1) = d;
    estDistance(i+1) = estd;
end

% Create error vectors for beacons of interest. i.e. Only the distance to
% the next beacon
% Times the robot reaches beacon: 5 13 21

beacons = [ 1 0;
            5 0;
            5 5;
            8 1];
figure(1)
plot(xpos_, ypos_);
hold on
plot(estimate(1,:),estimate(2,:));

plot(xpos, ypos)
plot(beacons(:,1),beacons(:,2), 'o')
%plot(x_tris, y_tris)
title("Robot path")
xlabel('X coordinate')
xlim([0 9])
ylabel('Y coordinate')
legend('Ideal path', 'EKF', 'triangulation', 'Beacons')
labels = {'Start','Beacon 1', 'Beacon 2', 'Beacon 3'};
text(beacons(:,1),beacons(:,2),labels,'VerticalAlignment','top','HorizontalAlignment','right')

figure(2)
plot(t,distance)
hold on
plot(t,estDistance)
title("Distance from Beacon Path")
xlabel('time (s)')
ylabel('Distance (m)')
legend('Noisy Raw Measurements', 'EKF Measurements')

