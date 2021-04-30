%%Other Estimation Method assignment
% We chose UKF (Unscented Kalman Filter)
clear all

%% System Simulation Assignment
b = .2; % Axle width (distance from left wheel to right wheel)
piv1 = .05 * pi / 3;
piv2 = -1 * (piv1 + piv1 * .927/1.57);
vl = [1 1 1 1 1 -piv1 -piv1 -piv1 1 1 1 1 1 -piv2 -piv2 -piv2 1 1 1 1 1 0 0]; % Measured velocities of left wheel
vr = [1 1 1 1 1  piv1  piv1  piv1 1 1 1 1 1  piv2  piv2  piv2 1 1 1 1 1 0 0]; % Measured velocities of right wheel
L = length(vl); % Number of samples in vl, vr (Currently 23)
t = 0:1:L; % Time vector 
% Summary of desired path:
% Start at origin
%  Arrive at beacon one (5,0) at k = 5
%  Pivot towards beacon two (5,5) for 3 seconds
%  Arrive at beacon two at k = 13
%  Pivot towards beacon three (8,1) for 3 seconds 
%  Arrive at beacon three at k = 21



%% Noise and disturbances section
mu = 0; % Mean for noise (unused)
sigma_xy = .15; % Covariance for measurement noise
sigma_xyDot = .05; % Covariance for process/plant noise
sigma_phi = .02; % Covariance for meas. noise (angle)
sigma_phiDot = .001; % Cov for process/plant noise

sigma_r = .05; % covar for Beacon distance meas.
sigma_th = .01; % covar for beacon angle meas.
sigma_phi2 = .001; % covar for vehic pose beacon meas.


% Measurement noise
x_noise = sigma_xy * randn(L+1,1); % add mu for non-zero mean
y_noise = sigma_xy * randn(L+1,1);
phi_noise = sigma_phi * randn(L+1,1);

d1_noise = sigma_r * randn(L+1,1); % noise for beacon distance 1
d2_noise = sigma_r * randn(L+1,1);
d3_noise = sigma_r * randn(L+1,1);
th1_noise = sigma_th * randn(L+1,1);
th2_noise = sigma_th * randn(L+1,1);
th3_noise = sigma_th * randn(L+1,1);
phi_noise2 = sigma_phi2 * randn(L+1,1);

% Process noise
xDot_noise = sigma_xyDot * randn(L+1,1);
yDot_noise = sigma_xyDot * randn(L+1,1);
phiDot_noise = sigma_phiDot * rand(L+1,1);
xDot_noise = zeros(L+1,1); %testing
yDot_noise = zeros(L+1,1);
phiDot_noise = zeros(L+1,1);


% Left disturbances
dist_l = zeros(L+1, 1);
%dist_l(5:6) = .2;

% Right disturbances
dist_r = zeros(L+1,1);
dist_r(17:20) = .15;

%% System initialization section
% Initial Parameters
Q = diag([sigma_xyDot sigma_xyDot sigma_phiDot]); % Covariance matrix for plant noise (?)
%R = diag([sigma_xy sigma_xy sigma_phi]); % Covariance matrix for measurement noise (?)
R = diag([sigma_r sigma_th sigma_r sigma_th sigma_r sigma_th sigma_phi2]);
P0 = eye(3,3)./500; % Covariance of initial state (~0 since known exactly)
x0 = [0;0;0]; % Initial coordinates

Pk_plus = P0;
xk_plus = x0;

x_tris = zeros(1,L+1);
y_tris = zeros(1,L+1);

p1 = [5; 0]; % beacon 1 location
p2 = [5; 5]; % beacon 2 location
p3 = [8; 1]; % beacon 3 location
% p_robot = zeros(L+1,2);

% Error measurement matrices
simError = zeros(L+1,3);
estError = zeros(L+1,3);

distance = zeros(L+1,1);
estDistance = zeros(L+1,1);

n = 3; % State Dimensions
m = 7; % Measurement Dimensions

e = eye(3,3);

avg = (1/(2*n)) * ones(2*n,1);

state_samples = zeros(n,2*n); % Sample distribution matrix
state_sim = zeros(n,2*n);   %This will store the simulated values state_sim = f(state_samples,u)
out_sim = zeros(m,2*n);   %This will store the simulated values out_sim = g(x_estimation_samples,u)

xDot=0;
yDot=0;
phiDot=0;

% Pos for "position", not "positive"
xpos_ =      zeros(1,L+1); % For ideal sim
xpos =       zeros(1,L+1); % After process noise and disturbances
xpos_noisy = zeros(1,L+1); % After measurement noise
ypos_ =      zeros(1,L+1); 
ypos =       zeros(1,L+1); 
ypos_noisy = zeros(1,L+1);
phipos =     zeros(1,L+1);
phipos_ =    zeros(1,L+1);
phipos_noisy = zeros(1,L+1);

estimate = zeros(3,L+1);

for i = 1:L
    %% Ideal sim section (no noise, no disturbance):
    xDot_ = .5 * cos(phipos_(i)) * (vr(i) + vl(i)); % Incremental change in x
    xpos_(i+1) = xDot_ + xpos_(i); % iterate sim
    
    yDot_ = .5 * sin(phipos_(i))* (vr(i) + vl(i));
    ypos_(i+1) = yDot_ + ypos_(i); % iterate simulation
    
    phiDot_ = vr(i)/b - vl(i)/b;
    phipos_(i+1) = phiDot_ + phipos_(i);
    
    %% Realistic sim section:
    xDot = .5 * cos(phipos(i)) * (vr(i) + dist_r(i) + vl(i) + dist_l(i)) + xDot_noise(i); %Includes disturbance and measurement noise
    xpos(i+1) = xDot + xpos(i); % Iterate simulation
    xpos_noisy(i+1) = xpos(i+1) + x_noise(i+1); % Collect measurement (with noise)

    yDot = .5 * sin(phipos(i))* (vr(i) + dist_r(i) + vl(i) + dist_l(i)) + yDot_noise(i);
    ypos(i+1) = yDot + ypos(i);
    ypos_noisy(i+1) = ypos(i+1) + y_noise(i+1);

    phiDot = 1/b * (vr(i) + dist_r(i) - (vl(i) + dist_l(i))) + phiDot_noise(i);
    phipos(i+1) = phiDot + phipos(i) ;
    phipos_noisy(i+1) = phipos(i+1) + phi_noise(i+1);

    % outk = [xpos_noisy(i+1); ypos_noisy(i+1); phipos_noisy(i+1)];

    % Vincent's requested output version:
    outk = [  sqrt( (p1(1)-xpos(i))^2 + (p1(2)-ypos(i))^2);
            atan2( p1(2)-ypos(i), p1(1) - xpos(i));
            
            sqrt( (p2(1)-xpos(i))^2 + (p2(2)-ypos(i))^2) ;
            atan2( p2(2)-ypos(i), p2(1) - xpos(i));
            
            sqrt( (p3(1)-xpos(i))^2 + (p3(2)-ypos(i))^2);
            atan2( p3(2)-ypos(i), p3(1) - xpos(i));
            
            phipos(i)];
        
     outk = outk + [d1_noise(i+1);
                    th1_noise(i+1);
                    d2_noise(i+1);
                    th2_noise(i+1);
                    d3_noise(i+1);
                    th3_noise(i+1);
                    phi_noise2(i+1);
                    ];

    %% UKF TIME UPDATE SECTION

    %     Get standard deviation matrix
    M = chol(Pk_plus);
    %     Simulate 2n points through state transition
    for j = 1:2*n
    %     Generate samples
        if(j<=n)
            state_samples(:,j) = xk_plus + sqrt(n)*M*e(:,j);
        else
            state_samples(:,j) = xk_plus - sqrt(n)*M*e(:,j-n);
        end
    %     Evaluate Simulated values
        state_sim(:,j) = xk_plus + [cos(state_samples(3,j))*vr(i)/2 + cos(state_samples(3,j))*vl(i)/2;
                                sin(state_samples(3,j))*vr(i)/2 + sin(state_samples(3,j))*vl(i)/2;
                                vr(i)/b - vl(i)/b];
    end

    % Take average mean as predicted mean
    x_kplus1_minus = state_sim*avg;

    for j = 1:2*n
        P_kplus1_minus = zeros(3,3);

    % Take average covariance as predicted covariance
        P_kplus1_minus = P_kplus1_minus + (1/(2*n)) * (state_sim(:,j)-x_kplus1_minus) * (state_sim(:,j)-x_kplus1_minus)';
    end
    P_kplus1_minus = P_kplus1_minus + Q;


    %% UKF MEASUREMENT UPDATE section

    xk_minus = x_kplus1_minus;
    Pk_minus = P_kplus1_minus;

    % Simulate 2n points through measurement
    for j = 1:2*n
    %     Generate samples
        if(j<=n)
            state_samples(:,j) = xk_minus + sqrt(n)*M*e(:,j);
        else
            state_samples(:,j) = xk_minus - sqrt(n)*M*e(:,j-n);
        end
    %     Evaluate Simulated output values (This is where Beacon measurements
    %     go)
         %out_sim(:,j) = state_samples(:,j);        % Just pass location

          % Get output from simulated beacon measurements
        out_sim(:,j) = [sqrt( (p1(1)-state_samples(1,j))^2 + (p1(2)-state_samples(2,j))^2);
                        atan2( p1(2)-state_samples(2,j), p1(1) - state_samples(1,j));
                        sqrt( (p2(1)-state_samples(1,j))^2 + (p2(2)-state_samples(2,j))^2) ;
                        atan2( p2(2)-state_samples(2,j), p2(1) - state_samples(1,j));
                        sqrt( (p3(1)-state_samples(1,j))^2 + (p3(2)-state_samples(2,j))^2);
                        atan2( p3(2)-state_samples(2,j), p3(1) - state_samples(1,j));
                        state_samples(3,j) ];
    end
    out_est = out_sim * avg;
    
    for j = 1:2*n
        P_xy = zeros(n,m);
        P_y = zeros(m,m);

        % Take average covariance as predicted covariance
        P_xy = P_xy + (1/(2*n)) * (state_samples(:,j)-xk_minus) * (out_sim(:,j)-out_est)';
        P_y = P_y + (1/(2*n)) * (out_sim(:,j)-out_est) * (out_sim(:,j)-out_est)';
    end
    P_y = P_y + R;

    % Calculate Kalman Gain
    K = P_xy / P_y;

    xk_plus = xk_minus + K*(outk - out_est);
    Pk_plus = Pk_minus - K*P_xy';
   
    %% Observation section
    estimate(:,i+1) = xk_plus;
    % Distance with noisy simulated measurements
    d_1 = norm([xpos_noisy(i+1);ypos_noisy(i+1)]-p1); % Distance to beacon 1
    d_2 = norm([xpos_noisy(i+1);ypos_noisy(i+1)]-p2); % Distance to beacon 2
    d_3 = norm([xpos_noisy(i+1);ypos_noisy(i+1)]-p3); % Distance to beacon 3
    simError(i+1, :) = [d_1, d_2, d_3];
        
    % Estimate error from actual simulated path (with noise)
    d = norm([xpos_noisy(i+1);ypos_noisy(i+1)]-[xpos(i+1);ypos(i+1)]); % raw measurement error to actual simulated path (with process noise and disturbances)
    estd = norm([xk_plus(1);xk_plus(2)]-[xpos(i+1);ypos(i+1)]); % Filter-estimated distance 
    distance(i+1) = d;
    estDistance(i+1) = estd;
end




%% Plotting section
beacons = [ 1 0;
            5 0;
            5 5;
            8 1];
figure(1)
plot(xpos_, ypos_); % Intended path (no process or measurement noise)
hold on
plot(xpos, ypos); % Actual path (no measurement noise)
plot(xpos_noisy, ypos_noisy) % Actual state (includes measurement and process noise)
plot(estimate(1,:),estimate(2,:)); % UKF results
plot(beacons(:,1),beacons(:,2), 'o')

title("Robot path")
xlabel('X coordinate')
xlim([0 9])
ylabel('Y coordinate')
legend('Intended Path', 'actual path', 'unfiltered measurement results', 'UKF results')
labels = {'Start','Beacon 1', 'Beacon 2', 'Beacon 3'};
text(beacons(:,1),beacons(:,2),labels,'VerticalAlignment','top','HorizontalAlignment','right')

figure(2)
plot(t,distance)
hold on
plot(t,estDistance)
title("Distance from Beacon Path")
xlabel('time (s)')
ylabel('Distance (m)')
legend('Noisy Raw Measurements', 'UKF Measurements')

