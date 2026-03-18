%======================================================================
% Sun‑centric Two‑Body Propagation and Diagnostics
%======================================================================
% This script demonstrates a simple Kepler‑only model for three bodies
% (Earth, the binary‑asteroid Didymos, and the DART impactor) orbiting the 
%  Sun.
%
%   Steps:
%   1.  Define heliocentric Cartesian state vectors (position [km] and velocity
%       [km/s]) for each body.
%   2.  Propagate the states using MATLAB's ODE45 solver and the
%       two‑body equations of motion (acceleration = -mu*r/|r|^3).
%       *   Didymos : 70 000 000 s (≈ 2 yr) in 100 000‑s steps
%       *   DART    : 26 000 000 s in 100 000‑s steps
%       *   Earth   : 70 000 000 s in 100 000‑s steps
%       Tight tolerances (RelTol/AbsTol = 1e‑13) keep numerical drift minimal.
%
%   3.  Plot the three 3‑D trajectories together for visual inspection.
%
%   4.  For Didymos only:
%        • Plot the magnitude of position, velocity, and gravitational
%          acceleration versus time.
%        • Compute and plot specific orbital energy (½v² – μ/r) – should be constant.
%        • Compute and plot specific angular‑momentum vector h = r × v and its
%          components – also conserved.
%        • Print the final Cartesian state after the integration.
%
%   5.  Convert every propagated Didymos state back to classical orbital
%      elements (e, a, i, ω, Ω, ν) using the routine
%      `cartesian_to_orbital`.  Plot each osculating element versus time;
%      for a pure two‑body problem the curves should be flat (aside from
%      tiny numerical noise).
%
%   6.  Helper functions:
%        • propagate_2BP – returns the time derivative [vx vy vz ax ay az] 
%          for the Sun‑centered two‑body problem.

%Inital Conditions
ic_earth = [
    6.82500E+07;   %x
    1.30864E+08;   %y
    1.81329E+04;   %z
   -2.67639E+01;   %xv
    1.38981E+01;   %yv
   -9.22794E-04    %zv
];

fc_earth = [
   -1.21652E+08;   
    8.56703E+07;   
   -2.05971E+03;   
   -1.74946E+01;   
   -2.42677E+01;   
   -3.76780E-03    
];

ic_didymos = [
   -2.39573E+08;   
   -2.35661E+08;  
    9.54384E+06;  
    1.24732E+01;   
   -9.74427E+00;  
   -8.78661E-01    
];

ic_dart = [
    6.82409E+07;
    1.30854E+08;   
    1.52197E+04;   
   -3.06997E+01;  
    8.11796E+00;   
    3.95772E+00   
];

%Propagate Didymos Orbit for 70,000,000 seconds
didymos_trange = [0:100000:70000000];
didymos_options = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);
[didymos_t, didymos_RV] = ode45(@propagate_2BP, didymos_trange, ic_didymos, didymos_options, 1.32712e11);

%Propagate DART oribt for 26,000,000 seconds
dart_trange = [0:100000:26000000];
dart_options = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);
[dart_t, dart_RV] = ode45(@propagate_2BP, dart_trange, ic_dart, dart_options, 1.32712e11);

% Propagate Earth Orbit for 70,000,000 seconds
earth_trange = [0:100000:70000000];
earth_options = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);
[earth_t, earth_RV] = ode45(@propagate_2BP, earth_trange, ic_earth, earth_options, 1.32712e11);

%Plot Didymos Orbit
figure;
hold on;
plot3(didymos_RV(:, 1), didymos_RV(:, 2), didymos_RV(:, 3), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Didymos Orbit');
plot3(ic_didymos(1), ic_didymos(2), ic_didymos(3), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b', 'DisplayName', 'Didymos Initial');


%Plot Dart orbit
plot3(dart_RV(:, 1), dart_RV(:, 2), dart_RV(:, 3), 'r-', 'LineWidth', 1.5, 'DisplayName', 'DART Orbit');
plot3(ic_dart(1), ic_dart(2), ic_dart(3), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'DART Initial');



%Plot Earth Oribt
plot3(earth_RV(:, 1), earth_RV(:, 2), earth_RV(:, 3), 'g-', 'LineWidth', 1.5, 'DisplayName', 'Earth Orbit');

%Axis and Legend
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
title('Earth, Didymos, and DART Orbits');
grid on;
axis equal;
legend('Location', 'Best');
hold off;


%Plot position, velocity, acceleration magnitudes of Didymos
figure;
subplot(3,1,1);
plot(didymos_t, vecnorm(didymos_RV(:, 1:3), 2, 2), 'b-', 'LineWidth', 1.5);
title('Didymos Orbit: Position, Velocity, and Acceleration Magnitudes');
ylabel('Position Magnitude (km)');
grid on;

subplot(3,1,2)
plot(didymos_t, vecnorm(didymos_RV(:, 4:6), 2, 2), 'g-', 'LineWidth', 1.5);
ylabel('Velocity Magnitude (km/s)');
grid on;

subplot(3,1,3)
acc = -1.32712e11 * didymos_RV(:, 1:3) ./ (vecnorm(didymos_RV(:, 1:3), 2, 2).^3);
acc_magnitude = vecnorm(acc, 2, 2);
plot(didymos_t, acc_magnitude, 'r-', 'LineWidth', 1.5);
ylabel('Acceleration Magnitude (km/s^2)');
xlabel('Time (seconds)');
grid on;

disp(['3. The magnitude plots make sense because given the position vector being squared in the denominator,' ...
    'as the value of r decreases, which happens as didymos gets closer to the sun, the acceleration will increase.'])

%Final State Didymos Oribt after 70,000,000 seconds
fs_didymos = didymos_RV(end,:);
fs_position = fs_didymos(1:3);
fs_velocity = fs_didymos(4:6);
disp('Final State of Didymos Orbit:');
disp('Position (km):');
disp(fs_position);
disp('Velocity (km/s):');
disp(fs_velocity);

% Calculate specific energy of didymos
speed = vecnorm(didymos_RV(:, 4:6), 2, 2);
dist = vecnorm(didymos_RV(:, 1:3), 2, 2);
energy = (0.5 * speed.^2) - 1.32712e11./dist;

%Plot specific ennergy of didymos over time
figure;
plot(didymos_t, energy, 'b', 'LineWidth', 1.5);
xlabel('Time (seconds)');
ylabel('Specific Energy (km^2/s^2)');
grid on;

% Calculate specific angular momentum of didymos
angular_momentum = cross(didymos_RV(:, 1:3), didymos_RV(:, 4:6));
angular_momentum_magnitude = vecnorm(angular_momentum, 2, 2);

%Plot specific angular momentum of didymos
figure;
subplot(4,1,1);
plot(didymos_t, angular_momentum_magnitude, 'g-', 'LineWidth', 1.5);
title('Specific Angular Momentum Magnitude of Didymos Orbit');
ylabel('Angular Momentum Magnitude (km^2/s)');
grid on;

subplot(4,1,2);
plot(didymos_t, angular_momentum(:, 1), 'r-', 'LineWidth', 1.5);
title('Specific Angular Momentum X Component');
ylabel('X Component (km^2/s)');
grid on;

subplot(4,1,3);
plot(didymos_t, angular_momentum(:, 2), 'b-', 'LineWidth', 1.5);
hold on;

title('Specific Angular Momentum Y Component');
ylabel('Y Component (km^2/s)');
xlabel('Time (seconds)');
legend('Y Component');
grid on;

subplot(4,1,4);
plot(didymos_t, angular_momentum(:, 3), 'm-', 'LineWidth', 1.5);
title('Specific Angular Momentum Z Component');
ylabel('Z Component (km^2/s)');
xlabel('Time (seconds)');
legend('Z Component');
grid on;

%Loop through each time step and compute orbital elements
num_steps = length(didymos_t);
elements = zeros(num_steps, 6);
for i = 1:num_steps
    r = didymos_RV(i, 1:3);
    v = didymos_RV(i, 4:6);
    [e, a, incl, w, omega, V] = cart_to_orb(r, v, 1.32712e11);
    elements(i, :) = [e, a, incl, w, omega, V];
end

% Plot osculating orbital elements using subplot
figure;

subplot(5,2,1);
plot(didymos_t, elements(:, 1), 'b-', 'LineWidth', 1.5);
title('Eccentricity (e)');
xlabel('Time (seconds)');
grid on;

subplot(5,2,2);
plot(didymos_t, elements(:, 2), 'r-', 'LineWidth', 1.5);
title('Semi-Major Axis (a)');
xlabel('Time (seconds)');
grid on;

subplot(5,2,3);
plot(didymos_t, elements(:, 3), 'g-', 'LineWidth', 1.5);
title('Inclination (i)');
xlabel('Time (seconds)');
grid on;

subplot(5,2,4);
plot(didymos_t, elements(:, 4), 'm-', 'LineWidth', 1.5);
title('Argument of Periapsis (w)');
xlabel('Time (seconds)');
grid on;

subplot(5,2,5);
plot(didymos_t, elements(:, 5), 'c-', 'LineWidth', 1.5);
title('Longitude of Ascending Node (\omega)');
xlabel('Time (seconds)');
grid on;

subplot(5,2,6);
plot(didymos_t, elements(:, 6), 'k-', 'LineWidth', 1.5);
title('True Anomaly (V)');
xlabel('Time (seconds)');
grid on;

function out = propagate_2BP(t,input, mu)
    x = input(1);
    y = input(2);
    z = input(3);
    vx = input(4);
    vy = input(5);
    vz = input(6);

    r = norm([x, y, z]);
    ax = -mu * x / r^3;
    ay = -mu * y / r^3;
    az = -mu * z / r^3;
    
    % Output vector
    out = [vx; vy; vz; ax; ay; az];
end