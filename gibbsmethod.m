function v_i = gibbsmethod(r1, r2, r3, mu, i)
% GIBBSMETHOD  Compute velocity at one of three coplanar position vectors.
%
%
%   INPUTS
%       r1, r2, r3 – 3‑vectors of position (km) at three successive epochs
%       mu – gravitational parameter (km^3/s^2)
%       i – index of the desired epoch (1, 2, or 3)
%
%   OUTPUT
%       v_i – 3‑vector velocity (km/s) at the selected epoch
%
%   The implementation follows the classic Gibbs method (see
%   Curtis, *Orbital Mechanics for Engineering Students*, Sec. 3.4).  It
%   assumes the three position vectors lie on the same orbit and are not
%   colinear.
%
%   The method forms three auxiliary vectors D, N and S, then evaluates
%   a scaling factor L = sqrt(mu/(|D| |N|)).  The final velocity is
%
%        v = (L / |r_i|) * B + L * S
%
%   where B = cross(D, r_i).


    r1_mag = norm(r1); 
    r2_mag = norm(r2);
    r3_mag = norm(r3);
    D = cross(r2,r3)+cross(r3,r1)+cross(r1,r2);
    D_mag = norm(D);
    N = r1_mag*(cross(r2,r3))+r2_mag*(cross(r3,r1))+r3_mag*(cross(r1,r2));
    N_mag = norm(N);
    S = (r1_mag-r2_mag)*r3+(r3_mag-r1_mag)*r2+(r2_mag-r3_mag)*r1;
 
 % Computation for velocity depending on i
    if i == 1
        B = cross(D,r1);
        L = sqrt(mu/(D_mag*N_mag));
        v_i = (L/r1_mag)*B+L*S;
    elseif i == 2
        B = cross(D,r2);
        L = sqrt(mu/(D_mag*N_mag));
        v_i = (L/r2_mag)*B+L*S;
    else
        B = cross(D,r3);
        L = sqrt(mu/(D_mag*N_mag));
        v_i = (L/r3_mag)*B+L*S;
    end

end