%======================================================================
% ESCAPE_THRUST_ANALYSIS
%======================================================================
%
% This script investigates how a spacecraft in a low‑Earth circular orbit
% (r0 = 8000 km) escapes when a constant, low‑thrust acceleration is
% applied along the +y‑axis.  Both an **analytical escape‑time formula**
% (derived from the energy equation for a constant thrust) and a 
% **numerical integration** of the two‑body equations with thrust are used.
%----------------------------------------------------------------------

%ICs
mu = 3.986*10^5;
T = [0 1e-4 0]; % kN/kg
r0 = [8000 0 0];

v0 = [0 0 0];
v0(2) = sqrt(mu/norm(r0));

aT = norm(T);
t_esc = (norm(v0)/aT)*(1-((20*aT^2*norm(r0)^2/norm(v0)^4)^(1/8)));

%ode45 options
tolerance = 1e-13;
options = odeset('RelTol', tolerance, 'AbsTol', tolerance);
input = [r0 v0];
tinitial = 0;
tfinal = 2*t_esc;
stepsize = 1;
trange = 0:stepsize:tfinal;

% Propagation
[t, RV] = ode45(@(t, R) propagate_2BP(R, mu, aT), trange, input, options);

% Plotting
rx = RV(:,1);
ry = RV(:,2);
figure(1)
plot(rx,ry)
axis equal
grid on
title('Trajectory in the xy‑plane (low‑thrust escape)')
xlabel('x  [km]')
ylabel('y  [km]')

% Final Velocity
vf = [RV(end,4) RV(end,5) RV(end,6)];

% Calculate escape velocity
r_esc = (norm(r0)*norm(v0))/(20*aT^2*norm(r0)^2)^(1/4);

%While loop to find esc_t
i = 1;
n = length(RV);
found = false;
escapetime = NaN;

while i <= n && ~found

    %Finds current r distance
    current_r = norm(RV(i,1:3));

    %If current r is greater than r_esc and store time
    if current_r > r_esc
        esc_t = t(i);
        found = true;
    end 
    i = i+1;
end

% Initial arrays and variables
step = 10;
trange_f = 0:step:10*t_esc;
aT_range = linspace(1e-5,1e-3,10);

%Initalize arrays to store ranging values 1-10
t_esc_a = zeros(1,10);
r_esc_a = zeros(1,10);
esc_t_a = zeros(1,10);

% Loop through each thrust value
for k = 1:length(aT_range)

    % Calculate escape radius and time for this thrust level
    r_esc_a(k) = (norm(r0) * norm(v0)) / (20 * aT_range(k)^2 * norm(r0)^2)^(1/4);
    t_esc_a(k) = (norm(v0) / aT_range(k)) * (1 - ((20 * aT_range(k)^2 * norm(r0)^2 / norm(v0)^4)^(1/8)));

    % Run the ODE solver for the current thrust value
    [t_thrust, RV_T] = ode45(@(t, R) propagate_2BP(R, mu, aT_range(k)), trange_f, input, options);

    % Initialize variables for while loop to find escape time
    i = 1;
    n = length(RV_T);
    found = false;
    esc_t_a(k) = NaN; 

    % While loop to find escape time
    while i <= n && ~found
            current_r = norm(RV_T(i, 1:3));

        if current_r > r_esc_a(k)
            esc_t_a(k) = t_thrust(i);
            found = true;
        end
    i = i + 1;

    end
end

% Plotting
figure(2)
loglog(aT_range,esc_t_a,'ro-','LineWidth',1.5,'MarkerSize',6); hold on
loglog(aT_range,t_esc_a,'b^--','LineWidth',1.5,'MarkerSize',6)
grid on
title('Escape time vs. constant thrust magnitude')
xlabel('Thrust acceleration a_T  [km/s^2]')
ylabel('Escape time  t_{esc}  [s]')
legend('Numerical (ode45)','Analytical','Location','best')
hold off

function rdot = propagate_2BP(R, mu, aT)


    r = R(1:3);
    v = R(4:6);

    %Gravitational accel (x,y,z)
    a_grav = -mu*r/norm(r)^3;

    %Vecloity unit vector
    v_hat = v/norm(v);

    %Thrust accel
    a_thrust = aT*v_hat;

    %Total accel
    a_tot = a_grav + a_thrust;

    %output
    rdot = [v; a_tot];
end