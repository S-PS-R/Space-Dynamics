function [e,a,i,w,Omega,nu] = cart_to_orb(r,v,mu)
% CART2ORB  Convert Cartesian state (r,v) to classical orbital elements.
%   [e,a,i,w,Omega,nu] = CART2ORB(r,v,mu) returns
%       e     – eccentricity
%       a     – semi‑major axis (km)      (NaN for parabolic)
%       i     – inclination (rad)
%       w     – argument of perigee (rad)
%       Omega – longitude of ascending node (rad)
%       nu    – true anomaly (rad)
%
%   Input vectors r and v must be 1×3 row vectors.
%   mu is the standard gravitational parameter (km^3/s^2).

    h = cross(r,v);
    h_norm = norm(h);
    e_vec = (1/mu)*cross(v,h) - r/norm(r);
    e = norm(e_vec);
    k = [0 0 1];
    n_vec = cross(k,h);
    n_norm = norm(n_vec);
    tol = 1e-12;


    if e < 1 - tol
        energy = norm(v)^2/2 - mu/norm(r);
        a = -mu/(2*energy);
    else
        a = NaN;
    end


    i = acos( h(3)/h_norm );


    if n_norm > tol
        Omega = acos( n_vec(1)/n_norm );
        if n_vec(2) < 0
            Omega = 2*pi - Omega;
        end
    else
        Omega = 0;
    end


    if e > tol && n_norm > tol 
        w = acos( dot(n_vec,e_vec)/(n_norm*e) );
        if e_vec(3) < 0
            w = 2*pi - w;
        end
    elseif e > tol && n_norm <= tol
        w = acos( e_vec(3)/e );
        if e_vec(2) < 0
            w = 2*pi - w;
        end
    else 
        w = 0; 
    end

    if e > tol
        nu = acos( dot(e_vec,r)/(e*norm(r)) );
        if dot(r,v) < 0
            nu = 2*pi - nu;
        end
    else
        if n_norm > tol
            nu = acos( dot(n_vec,r)/(n_norm*norm(r)) );
            if dot(r,v) < 0
                nu = 2*pi - nu;
            end
        else
            nu = acos( r(1)/norm(r) );
            if r(2) < 0
                nu = 2*pi - nu;
            end
        end
    end
end