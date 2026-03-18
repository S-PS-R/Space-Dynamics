% ================================================================
% Kepler Prediction Problem – universal‑variable propagator
% ------------------------------------------------
% All quantities are in km, km/s, seconds.
% ------------------------------------------------
% Input: orbital elements (a, e, i, Ω, ω, ν0) and gravitational constant μ.
% Output: propagated state (r, v) after a user‑specified time T.
% ================================================================

%%{
a = 15000;
v0 = 145;
e = 0.4;
i = 60;
omega = 45;
w = 0;
u = 3.986*10^5;
%%}


%{
a = 21000;
v0 = 35;
e = 0.4;
i = 90;
omega = 0;
w = 0;
u = 3.986*10^5;
%}

v0 = deg2rad(v0);
i = deg2rad(i);
w = deg2rad(w);
omega = deg2rad(omega);


%Initial Radius 
[r0_v, V0_v] = orb_to_cart(u,e,a,i,w,omega,v0);
P = a*(1-e^2);
r0 = P/(1+e*cos(v0));
V0 = norm(V0_v);


%Algorithm
T = 3600;
Xn = sqrt(u)/a * T;
tol = 1*10^-6;
t = 5000;
X = Xn;

while abs(t-T) > 10

    psi = X^2/a;

    if psi > tol
        C = (1 - cos(sqrt(psi)))/psi;
        S = (sqrt(psi) - sin(sqrt(psi)))/sqrt(psi^3);
    
    elseif psi < -tol
        C = (1 - cosh(sqrt(-psi)))/psi;
        S = (sinh(sqrt(-psi))-sqrt(-psi))/sqrt((-psi)^3);

    else
        C = 1/2;
        S = 1/6;
       
    end

    r = X^2*C + (dot(r0_v,V0_v))/sqrt(u)  *  X*(1-psi*S) +   r0*(1-psi*C);
    dxdt = sqrt(u)/r;

    t = 1/sqrt(u) * (X^3)*S + (dot(r0_v,V0_v))/sqrt(u)*X^2*C + r0*X*(1-psi*S)
    X = X + (T-t)*dxdt

end 

f = 1 - X^2/r0 * C;
f_dot = sqrt(u)/(r0*r)*X*(psi*S-1);
g = t- X^3/sqrt(u) * S;
g_dot = 1- X^2/r * C;

r = f*r0_v + g*V0_v
v = f_dot*r0_v + g_dot*V0_v
