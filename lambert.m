
function [v0, vf, rp, e] = lambert(r0, rf, TOF, DM, mu)
% LAMBERT  Solve Lambert’s problem using the universal‑variable formulation.
%
%   INPUTS
%       r0   – 3‑vector, initial position (km)
%       rf   – 3‑vector, final   position (km)
%       TOF  – time‑of‑flight (s)
%       DM   – direction of motion (+1 for prograde, –1 for retrograde)
%       mu   – gravitational parameter (km^3/s^2)
%
%   OUTPUTS
%       v0   – 3‑vector, departure velocity (km/s)
%       vf   – 3‑vector, arrival   velocity (km/s)
%       rp   – periapsis radius of the transfer orbit (km)
%       e    – eccentricity of the transfer orbit (unitless)
%
%   The algorithm follows the universal‑variable method described in
%   Vallado, *Fundamentals of Astrodynamics and Applications* (4th ed.).
%   It iterates on the universal anomaly ψ until the computed time‑of‑flight
%   matches the requested TOF within a tolerance of 1 µs.
%
%   NOTE: The function now guards against division by zero when the
%   Stumpff functions C2 or C3 become singular (ψ → 0).
%----------------------------------------------------------------------------

    % Calculate the change in true anomaly
    r0_norm = norm(r0);
    rf_norm = norm(rf);
    delta_nu = acos(dot(r0, rf) / (r0_norm * rf_norm));
    A = DM * sqrt(r0_norm * rf_norm * (1 + cos(delta_nu)));
    if abs(A) < 1e-12
        error('Lambert:DegenerateCase', ...
              'Parameter A is zero – the problem is ill‑posed.');
    end
    
    % Inital guesses for universal anomaly and Stumpff functions
    psi = 0.0;
    C2 = 0.5;
    C3 = 1.0/6.0;
    psi_up = 4*pi^2;
    psi_low = -4*pi^2;
    
    tol = 1e-6;
    t_guess = 0.0;
    
    while abs(TOF - t_guess) > tol
        % Compute y(ψ)
        y = r0_norm + rf_norm + A * (psi*C3 - 1) / sqrt(C2);
        if y < 0   % y must be positive – adjust lower bound
            psi_low = psi + 0.1*(psi_up - psi);
        end
    
        % Compute x = sqrt(y/C2) safely
        if C2 <= 0
            % For ψ ≈ 0 we fallback to the series expansion of C2
            C2 = 0.5;
        end
        x = sqrt(y / C2);
    
        % Time‑of‑flight for current ψ
        t_guess = (x^3) * C3 + A * sqrt(y);
        t_guess = t_guess / sqrt(mu);
    
        % Update bounds
        if t_guess <= TOF
            psi_low = psi;
        else
            psi_up = psi;
        end
    
        % New trial ψ (bisection)
        psi = (psi_up + psi_low) / 2;
    
        % Update Stumpff functions for the new ψ
        if psi >  tol               % positive ψ  → elliptic
            sqrt_psi = sqrt(psi);
            C2 = (1 - cos(sqrt_psi)) / psi;
            C3 = (sqrt_psi - sin(sqrt_psi)) / (psi*sqrt_psi);
        elseif psi < -tol           % negative ψ → hyperbolic
            sqrt_mpsi = sqrt(-psi);
            C2 = (cosh(sqrt_mpsi) - 1) / (-psi);
            C3 = (sinh(sqrt_mpsi) - sqrt_mpsi) / ((-psi)*sqrt_mpsi);
        else                        % ψ ≈ 0 → use series expansion
            C2 = 0.5;
            C3 = 1.0/6.0;
        end
    end
    
    %Calculating lagarange coefficicents
    f = 1 - (y/norm(r0));
    g = A * sqrt(y/mu);
    g_dot = 1 - y / norm(rf);
    
    v0 = (rf - f*r0) /g;
    vf = (g_dot*rf - r0) /g;
    
    h = cross(r0, v0);
    e_vector = cross(v0, h)/mu -r0/norm(r0);
    e = norm(e_vector);
    
    epsilon = norm(v0)^2/2 -mu/norm(r0);
    
    a = -mu / (2 * epsilon);
    rp = a * (1 - e);

end