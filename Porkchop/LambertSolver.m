function [vi,vf]=LambertSolver(ro,rf,tof,tm,mu)
rom = norm(ro);
rfm = norm(rf);
tol = 1e-6;

dnu = acos(dot(ro,rf)/(rom*rfm));
A = tm*sqrt(rom*rfm*(1+cos(dnu)));
psi = 0;
c2 = 1/2;
c3 = 1/6;
psi_up = 4*pi^2;
psi_low = -4*pi^2;
tn = 0;
count = 0;
while abs(tn-tof) > tol
    yn = rom+rfm+(A*(psi*c3-1))/sqrt(c2);
    if A > 0 & yn < 0
        while yn < 0
            psi_low = psi_low+pi/4;
            psi = (psi_low+psi_up)/2;
            if psi > tol
                c2 = (1-cos(sqrt(psi)))/psi;
                c3 = (sqrt(psi)-sin(sqrt(psi)))/sqrt(psi^3);
            elseif psi < -tol
                c2 = (1-cosh(sqrt(-psi)))/psi;
                c3 = (sinh(sqrt(-psi))-sqrt(-psi))/sqrt(-psi^3);
            else
                c2 = 1/2;
                c3 = 1/6;
            end
            yn = rom+rfm+(A*(psi*c3-1))/sqrt(c2);
        end
    end
    if count >= 250 %if stuck in the while loop then break
        break
    end
    xn = sqrt(yn/c2);
    tn = (xn^3*c3+A*sqrt(yn))/sqrt(mu);
    
    if tn <= tof
        psi_low = psi;
    else
        psi_up = psi;
    end
    zn1 = (psi_low+psi_up)/2;
    
    if zn1 > tol
        c2 = (1-cos(sqrt(zn1)))/zn1;
        c3 = (sqrt(zn1)-sin(sqrt(zn1)))/sqrt(zn1^3);
    elseif zn1 < -tol
        c2 = (1-cosh(sqrt(-zn1)))/zn1;
        c3 = (sinh(sqrt(-zn1))-sqrt(-zn1))/sqrt(-zn1^3);
    else
        c2 = 1/2;
        c3 = 1/6;
    end
    psi = zn1;
    count = count+1;
end


f = 1-yn/rom;
g = A*sqrt(yn/mu);
gdot = 1-yn/rfm;
vi = (rf-f*ro)/g;
vf = (gdot*rf-ro)/g;