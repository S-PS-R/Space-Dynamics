function [r,v] = orb_to_cart(u,e,a,i,w,omega,V)
% orb_to_cart  Convert Keplerian orbital elements to Cartesian state vectors.
%   [r,v] = ORB_TO_CART(mu,e,a,i,w,Omega,nu) returns the position vector **r**
%   (km) and velocity vector **v** (km/s) in the inertial (ECI) frame for a
%   two‑body orbit described by the classical orbital elements:
%
%       mu    – gravitational parameter (km^3/s^2)   (GM of the central body)
%       e     – eccentricity (dimension‑less)
%       a     – semi‑major axis (km)
%       i     – inclination (rad)
%       w     – argument of perigee (rad)
%       Omega – longitude of the ascending node (rad)
%       nu    – true anomaly (rad)
%
%   Input vectors **r** and **v** are returned as **column vectors** (3×1)

    R3_w = [cos(-w) sin(-w) 0;-sin(-w)  cos(-w) 0;0 0   1];
    R1_i = [1   0   0;0 cos(-i) sin(-i);0   -sin(-i)    cos(-i)];
    R3_omega = [cos(-omega) sin(-omega) 0;-sin(-omega)  cos(-omega) 0;0 0   1];
    
    r_pqw = a*(1-e^2)/(1+e*cos(V));
    r_pqw_vector = [r_pqw*cos(V);    r_pqw*sin(V);    0];
    
    P = a*(1-e^2);
    v_pqw = sqrt(u/P)*[-sin(V);  e+cos(V);    0];
    
    T = R3_omega*R1_i*R3_w;
    
    r = T*r_pqw_vector;
    v = T*v_pqw;

end