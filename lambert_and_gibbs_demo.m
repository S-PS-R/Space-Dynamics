%======================================================================
% LAMBERT_AND_GIBBS_DEMO
%======================================================================
%
% This script demonstrates three classic orbital‑mechanics tools:
%
%  1) Lambert’s problem – compute the transfer velocities that connect
%     two position vectors in a given time‑of‑flight.
%  2) Numerical propagation of the Lambert solution with a pure
%     two‑body model (no thrust) using ode45, and a comparison with the
%     analytical arrival state.
%  3) Gibbs’ method – an initial‑orbit‑determination technique that
%     determines the velocity at the middle observation from three
%     position vectors.
%---------------------------------------------------------------------

%Case 1
r0 = [8000 0 0];
rf = [7000 7000 0];
TOF  = 3600;
DM = 1;
mu = 3.986*10^5;
disp('Case 1:')
[v0,vf,rp1,e1] = lambert(r0,rf,TOF,DM,mu);

% Case 2
r1 = [0.5 0.6 0.7]*6371;
r2 = [0 -1 0]*6371;
TOF2 = 16135;
DM2 = -1;
disp('Case 2:')
[v1,v2,rp2,e2] = lambert(r1,r2,TOF2,DM2,mu);

%2 body prop of case 2
r = r1;
v = v1;
input = [r v];

tolerance = 1E-13;
options = odeset('RelTol', tolerance, 'AbsTol', tolerance);
trange = 0:1:TOF2;
 
[~, RV] = ode45(@(t, R) propagate_2BP(t, R, mu), trange, input, options);

rx = RV(:, 1);
ry = RV(:, 2);
rz = RV(:, 3);
position = [rx(end), ry(end), rz(end)];

vx = RV(:, 4);
vy = RV(:, 5);
vz = RV(:, 6);
velocity = [vx(end), vy(end), vz(end)];


diff_r = r2 - position;
disp("Difference in final position")
disp(diff_r)

diff_v = v2 - velocity; 
disp("Difference in final velocity")
disp(diff_v)

figure;
plot3(rx, ry, rz, 'b');
hold on;
plot3(input(1), input(2), input(3), 'ro');
text(input(1), input(2), input(3), ' Initial');
plot3(rx(end), ry(end), rz(end), 'go');
text(rx(end), ry(end), rz(end), ' Final');
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Transfer Trajectory');
grid on;
axis equal;
hold off;

% Gibbs Method
clear; clc;
r1 = [-6.97949190E+03 2.08846535E+02 -5.73801140E+02];
r2 = [-6.96633930E+03 2.67474976E+02 -7.34881458E+02];
r3 = [-6.94996267E+03 3.25979638E+02 -8.95621695E+02];

mu = 3.986*10^5;
i=2;

disp('Velocity at Epoch 2')
V2  = gibbsmethod(r1,r2,r3,mu,i)
disp('Velocity Magnitude at Epoch 2')
disp(norm(V2))

function out = propagate_2BP(~, input, mu)
    % input = [rx ry rz vx vy vz]
    r = input(1:3);
    v = input(4:6);
    rdot = v;
    vdot = -mu / norm(r)^3 * r;
    out = [rdot; vdot];
end