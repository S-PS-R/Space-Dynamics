%======================================================================
% RIGID_BODY_ROTATION_AND_QUATERNION_DEMO
%======================================================================
% This script is a compact demonstration that ties together three
% (but related) topics in classical mechanics:
%
%   1. **Rotation‑matrix → axis‑angle → quaternion conversion**  
%   2. **Torque‑free rigid‑body dynamics (Euler’s equations)**  
%   3. **Geometric visualisation**  


R3 = [cosd(30) sind(30) 0; -sind(30) cosd(30) 0; 0 0 1];
R2 = [cosd(40) 0 -sind(40);0 1 0; sind(40) 0 cosd(40)];
R1 = [1 0 0; 0 cosd(10) sind(10); 0 -sind(10) cosd(10)];

R3_ = [cosd(-30) sind(-30) 0; -sind(-30) cosd(-30) 0; 0 0 1];
R2_ = [cosd(-40) 0 -sind(-40);0 1 0; sind(-40) 0 cosd(-40)];
R1_ = [1 0 0; 0 cosd(-10) sind(-10); 0 -sind(-10) cosd(-10)];

C1 = R3_*R2_*R1_


C = R1*R2*R3
phi = acosd(0.5*(C(1,1)+C(2,2)+C(3,3)-1))
e = 1/(2*sind(phi))*[C(2,3)-C(3,2);C(3,1)-C(1,3);C(1,2)-C(2,1)]

q1 = e(1)*sind(phi/2)
q2 = e(2)*sind(phi/2)
q3 = e(3)*sind(phi/2)
q4 = cosd(phi/2)

B0 = q4;
B1 = q1;
B2 = q2;
B3 = q3;
w = [0; 0.1; 0.2; 0]
B_dot = 0.5*[B0 -B1 -B2 -B3; B1 B0 -B3 B2; B2 B3 B0 -B1; B3 -B2 B1 B0]*w

I = [10 0 0;0 20 0;0 0 30];

%2a
w = [10 0 30]'*pi/180;

H = norm(I*w)
T = 0.5*transpose(w)*I*w

%2b
trange = [0 100];
tol = 1e-13;

options = odeset('RelTol',tol);
[t,w] = ode45(@angle,trange,w,options);

figure(1);
hold on;
plot(t,w(:,1));
plot(t,w(:,2));
plot(t,w(:,3));
legend('\omega1','\omega2','\omega3')
xlabel('Time (s)')
ylabel('Angular Velocity (rad/s)')


%Loop to calculate magnitude at time step
for i = 1:685

    %Calculate H at every time step
    H_mag(i) = norm(I*w(i,:)');
    T_mag(i) = .5*w(i,:)*I*w(i,:)';

    %Calculate difference 
    H_dev(i) = H_mag(i) - H;
    T_dev(i) = T_mag(i) - T;

    %Calculate time varying H for part d
    L = I*w(i,:)';
    H_vary(1,i) = L(1);
    H_vary(2,i) = L(2);
    H_vary(3,i) = L(3);
    
end


%Plot 2c
figure(2);
plot(t,H_dev)
xlabel('Time (s)');
ylabel('Angular momentum Deviation (kg*m^2/s')
title('Angular momentum deviation vs time')

figure(3)
plot(t,T_dev)
xlabel('Time (s)');
ylabel('Kinetic energy deviation (kg*m^2/s^2')
title('Kinetic energy deviation vs time')

%2d

%Sphere - Angular momentum
[H1 H2 H3] = sphere;

%Ellipsoid - Kinetic Energy
H4 = sqrt(2*I(1,1)*T);
H5 = sqrt(2*I(2,2)*T);
H6 = sqrt(2*I(3,3)*T);
[x,y,z] = ellipsoid(0,0,0,H4,H5,H6);

figure(4);
hold on
surf(H1*H,H2*H,H3*H,'FaceColor','c')
surf(x,y,z,'FaceColor','b');
plot3(H_vary(1,:),H_vary(2,:),H_vary(3,:),'.','Color','r')
legend('Sphere (H)', 'Ellipsoid (KE)', 'Time Varying H')
axis equal
xlabel('H1')
ylabel('H2')
zlabel('H3')
view([20,30])
hold off

clear

%2e
I = [10 0 0;0 20 0;0 0 30];
w0 = [1 15 0]'*pi/180;

trange = [0 100];
tol = 1e-13;
options = odeset('RelTol',tol);
[t,w] = ode45(@angle,trange,w0,options);

figure(5);
hold on;
plot(t,w(:,1));
plot(t,w(:,2));
plot(t,w(:,3));
legend('\omega1','\omega2','\omega3')
xlabel('Time (s)')
ylabel('Angular Velocity (rad/s)')

%2f
H0 = norm(I*w0);
T0 = 0.5*transpose(w0)*I*w0;

%Sphere - Angular momentum
[H1, H2, H3] = sphere;

%Ellipsoid - Kinetic Energy
H4 = sqrt(2*I(1,1)*T0);
H5 = sqrt(2*I(2,2)*T0);
H6 = sqrt(2*I(3,3)*T0);
[x,y,z] = ellipsoid(0,0,0,H4,H5,H6);

for i = 1:169

    %Calculate time varying H
    L = I*w(i,:)';
    H_vary(1,i) = L(1);
    H_vary(2,i) = L(2);
    H_vary(3,i) = L(3);
   
end


figure(6);
hold on
surf(H1*H0,H2*H0,H3*H0,'FaceColor','c')
surf(x,y,z,'FaceColor','b');
plot3(H_vary(1,:),H_vary(2,:),H_vary(3,:),'.','Color','r')
legend('Sphere (H)', 'Ellipsoid (KE)', 'Time Varying H')
axis equal
xlabel('H1')
ylabel('H2')
zlabel('H3')
view([20,40])
hold off


%2b Function
function [w_dot] = angle(t,w)

I = [10 0 0;0 20 0;0 0 30];

w1 = w(1);
w2 = w(2);
w3 = w(3);

%Calculate w_dot for ode45
w_dot(1,1) = -(I(3,3)-I(2,2))*w2*w3/I(1,1);
w_dot(2,1) = -(I(1,1)-I(3,3))*w3*w1/I(2,2);
w_dot(3,1) = -(I(2,2)-I(1,1))*w1*w2/I(3,3);

end

