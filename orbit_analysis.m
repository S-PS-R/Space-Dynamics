% -------------------------------------------------------------------------
% orbit_analysis.m
% -------------------------------------------------------------------------
% This script visualises a Keplerian orbit and shows how the osculating
% orbital elements vary when the state is repeatedly converted back and
% forth between Cartesian and orbital‑element representations.
%
% Usage:
%   1. Make sure orb_to_cart.m and cart_to_orb.m are on the MATLAB path.
%   2. Run the script:
%          orbit_demo          % or click the file in the Editor and press F5
%
%   The script will generate three figures:
%      • 3‑D trajectory of the orbit
%      • Energy‑deviation plot (should be zero for an ideal two‑body case)
%      • Six sub‑plots of the osculating orbital elements vs. time
% -------------------------------------------------------------------------

a = 20000;
e = 0.4;
i = 100*pi/180;
omega = 30*pi/180;
w = 15*pi/180;
u = 3.986*10^5;

period = 2*pi*sqrt(a^3/u);

figure;
num_points = 1000;
orbit_points = zeros(3, num_points);

for k = 1:num_points
    V = 2*pi*(k-1)/(num_points-1);
    [r, v] = orb_to_cart(u, e, a, i, w, omega, V);
    orbit_points(:, k) = r;
end

plot3(orbit_points(1, :), orbit_points(2, :), orbit_points(3, :), 'LineWidth', 1.5);
grid on;
title('Orbit in 3D');
xlabel('X');
ylabel('Y');
zlabel('Z');
axis equal;

[r0,v0] = orb_to_cart(u, e, a, i, w, omega, 15*pi/180);
[rf,vf] = orb_to_cart(u, e, a, i, w, omega, (360-15-1)*pi/180);

% Plot the deviation of energy
E0 = norm(v0)^2/2 - u/norm(r0);
Ef = norm(vf)^2/2 - u/norm(rf);
tspan = 0:1:1000;
figure;
hold on;
plot(tspan,(Ef-E0)*ones(size(tspan)),'LineWidth',1.5);
grid on;
xlabel('Time')
ylabel("Specific Energy")
title('Energy Deviation')

%Plot osculating orbital elements
semimajor_axis = zeros(size(V));
eccentricity = zeros(size(V));
inclination = zeros(size(V));
argument_of_periapsis = zeros(size(V));
longitude_of_ascending_node = zeros(size(V));
true_anomaly = zeros(size(V));

V = linspace(0, 2*pi*2, num_points);
for j = 1:length(V)
    [r, v] = orb_to_cart(u, e, a, i, w, omega, V(j));
    [eo, ao, io, wo, omegao, Vo] = cart_to_orb(r, v, u);
    semimajor_axis(j) = ao;
    eccentricity(j) = eo;
    inclination(j) = io;
    argument_of_periapsis(j) = wo;
    longitude_of_ascending_node(j) = omegao;
    true_anomaly(j) = Vo;
end

% Plot osculating orbital elements
figure;

subplot(2, 3, 1);
plot(V*500, semimajor_axis, 'LineWidth', 1.5);
grid on;
xlabel('Time');
ylabel('Semimajor Axis');
title('Osculating Semimajor Axis');

subplot(2, 3, 2);
plot(V*500, eccentricity, 'LineWidth', 1.5);
grid on;
xlabel('Time');
ylabel('Eccentricity');
title('Osculating Eccentricity');

subplot(2, 3, 3);
plot(V*500, inclination, 'LineWidth', 1.5);
grid on;
xlabel('Time');
ylabel('Inclination');
title('Osculating Inclination');

subplot(2, 3, 4);
plot(V*500, argument_of_periapsis, 'LineWidth', 1.5);
grid on;
xlabel('Time');
ylabel('Argument of Periapsis');
title('Osculating Argument of Periapsis');

subplot(2, 3, 5);
plot(V*500, longitude_of_ascending_node, 'LineWidth', 1.5);
grid on;
xlabel('Time');
ylabel('Longitude of Ascending Node');
title('Osculating Longitude of Ascending Node');

subplot(2, 3, 6);
plot(V*500, true_anomaly, 'LineWidth', 1.5);
grid on;
xlabel('Time');
ylabel('True Anomaly');
title('Osculating True Anomaly');