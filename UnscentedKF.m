% Here will be a simulation for the Unscented Kalman Filter
clear all

%% System Simulation Assignment
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



L = 23; % Number of samples in vl, vr
mu = 0; % Mean for noise in x and y measurements
sigma_xy = .15; % Covariance (or standard deviation?)
sigma_xyDot = .05;
sigma_phiDot = .01;
sigma_phi = .05; % Covariance for angle noise
x_noise = sigma_xy * randn(L,1); % add mu for non-zero mean
y_noise = sigma_xy * randn(L,1);
phi_noise = sigma_phi * randn(L,1);

xDot_noise = sigma_xyDot * randn(L,1);
yDot_noise = sigma_xyDot * randn(L,1);
phiDot_noise = sigma_phiDot * rand(L,1);


% Left disturbances
dist_l = zeros(23, 1);
dist_l(5:6) = .2;

% Right disturbances
dist_r = zeros(23,1);
dist_r(17:20) = .15;

% Initial Parameters
Q = diag([sigma_xyDot sigma_xyDot sigma_phiDot]); % Covariance matrix for plant noise (?)
R = diag([sigma_xy sigma_xy sigma_phi]); % Covariance matrix for measurement noise (?)
P0 = eye(3,3)./500; % Covariance of initial state (0 since known exactly)
x0 = [0;0;0]; % Initial coordinates

Pk_plus = P0;
xk_plus = x0;




p1 = [5; 0]; % beacon 1 location
p2 = [5; 5]; % beacon 2 location
p3 = [8; 1]; % beacon 3 location
p_robot = zeros(23,2);

n = 3; % State Dimensions
m = 3; % Measurement Dimensions

e = eye(3,3);

avg = (1/(2*n)) * ones(2*n,1);

x_samples = zeros(n,2*n); %Sample distribution matrix
x_sim = zeros(n,2*n);   %This will store the simulated values x_sim = f(x_samples,u)
y_sim = zeros(m,2*n);   %This will store the simulated values y_sim = g(x_estimation_samples,u)

xDot=0;
yDot=0;
phiDot=0;

xpos_noisy = zeros(1,23);
xpos = xpos_noisy;
ypos_noisy = zeros(1,23);
ypos = ypos_noisy;
phiPos = zeros(1,23);

estimate = zeros(3,23);

xlast = 0;
ylast = 0;
phiLast = 0;

for i = 0:L-1
    
%     Generate "actual" measurement
xDot = cos(phiLast)*vr(i+1)/2 + cos(phiLast)*vl(i+1)/2;
xlast = xDot + xlast ;
xpos(i+1) = xlast;
xpos_noisy(i+1) = xlast+ x_noise(i+1);

yDot = sin(phiLast)*vr(i+1)/2 + sin(phiLast)*vl(i+1)/2;
ylast = yDot + ylast;
ypos(i+1) = ylast ;
ypos_noisy(i+1) = ylast + y_noise(i+1);

phiDot = vr(i+1)/b - vl(i+1)/b;
phiLast = phiDot + phiLast ;
phiPos(i+1) = phiLast+ phi_noise(i+1);

yk = [xpos_noisy(i+1);ypos_noisy(i+1);phiPos(i+1)];



    
% TIME UPDATE

%     Get standard deviation matrix
M = chol(Pk_plus);
%     Simulate 2n points through state transition
for j = 1:2*n
%     Generate samples
    if(j<=n)
        x_samples(:,j) = xk_plus + sqrt(n)*M*e(:,j);
    else
        x_samples(:,j) = xk_plus - sqrt(n)*M*e(:,j-n);
    end
%     Evaluate Simulated values
    x_sim(:,j) = xk_plus + [cos(x_samples(3,j))*vr(i+1)/2 + cos(x_samples(3,j))*vl(i+1)/2;
                            sin(x_samples(3,j))*vr(i+1)/2 + sin(x_samples(3,j))*vl(i+1)/2;
                            vr(i+1)/b - vl(i+1)/b];
end

% Take average mean as predicted mean
x_kplus1_minus = x_sim*avg;

for j = 1:2*n
    P_kplus1_minus = zeros(3,3);

% Take average covariance as predicted covariance
    P_kplus1_minus = P_kplus1_minus + (1/(2*n)) * (x_sim(:,j)-x_kplus1_minus) * (x_sim(:,j)-x_kplus1_minus)';
end
P_kplus1_minus = P_kplus1_minus + Q;


% MEASUREMENT UPDATE

xk_minus = x_kplus1_minus;
Pk_minus = P_kplus1_minus;

%     Simulate 2n points through measurement
for j = 1:2*n
%     Generate samples
    if(j<=n)
        x_samples(:,j) = xk_minus + sqrt(n)*M*e(:,j);
    else
        x_samples(:,j) = xk_minus - sqrt(n)*M*e(:,j-n);
    end
%     Evaluate Simulated output values (This is where Beacon measurements
%     go)
    y_sim(:,j) = x_samples(:,j);        %Just pass location
end
    y_est = y_sim * avg;
    
for j = 1:2*n
    P_xy = zeros(n,m);
    P_y = zeros(m,m);

% Take average covariance as predicted covariance
    P_xy = P_xy + (1/(2*n)) * (x_samples(:,j)-xk_minus) * (y_sim(:,j)-y_est)';
    P_y = P_y + (1/(2*n)) * (y_sim(:,j)-y_est) * (y_sim(:,j)-y_est)';
end
P_y = P_y + R;

% Calculate Kalman Gain
K = P_xy /(P_y);

xk_plus = xk_minus + K*(yk - y_est);
Pk_plus = Pk_minus - K*P_xy';


 estimate(:,i+1) = xk_plus;
%     Distance with noisy simulated measurements
    d_1 = norm([xpos_noisy(i+1);ypos_noisy(i+1)]-p1); % Distance to beacon 1
    d_2 = norm([xpos_noisy(i+1);ypos_noisy(i+1)]-p2); % Distance to beacon 2
    d_3 = norm([xpos_noisy(i+1);ypos_noisy(i+1)]-p3); % Distance to beacon 3
    simError(i+1, :) = [d_1, d_2, d_3];
    
    
%     Distance with EKF estimations
    estd_1 = norm([xk_plus(1);xk_plus(2)]-p1); % Distance to beacon 1
    estd_2 = norm([xk_plus(1);xk_plus(2)]-p2); % Distance to beacon 2
    estd_3 = norm([xk_plus(1);xk_plus(2)]-p3); % Distance to beacon 3
    estError(i+1, :) = [estd_1 estd_2 estd_3];
    
    
    
%     Distance from path
    d = norm([xpos_noisy(i+1);ypos_noisy(i+1)]-[xpos(i+1);ypos(i+1)]); % Distance to beacon path
    estd = norm([xk_plus(1);xk_plus(2)]-[xpos(i+1);ypos(i+1)]); % Distance to beacon path
    distance(i+1) = d;
    estDistance(i+1) = estd;
    
    
%     Triangulation
% alpha = atan(1/3);
% gamma = acos( (estd_1^2 +10 - estd_2^2) / (2* estd_1 * sqrt(10)));
% 
% x_robot = p1(1) + estd_1*cos(gamma+alpha);
% y_robot = p1(2) + estd_1*sin(gamma+alpha);
% 
% p_robot(i+1, :) = [x_robot y_robot];

end
% Create error vectors for beacons of interest. i.e. Only the distance to
% the next beacon
% Times the robot reaches beacon: 5 13 21

% simDistance = [simError(1:5,1);
%                 simError(6:13,2);
%                 simError(14:23,3)];
% estDistance = [estError(1:5,1);
%                 estError(6:13,2);
%                 estError(14:23,3)];



beacons = [ 1 0;
            5 0;
            5 5;
            8 1];
figure(1)
plot(xpos,ypos);
hold on
plot(estimate(1,:),estimate(2,:));

plot(xpos_noisy, ypos_noisy)
plot(beacons(:,1),beacons(:,2), 'o')
title("Robot path")
xlabel('X coordinate')
xlim([0 9])
ylabel('Y coordinate')
legend('raw results', 'UKF results', 'noisy results')
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

