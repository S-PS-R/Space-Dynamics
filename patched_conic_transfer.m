%======================================================================
% Mercury‑to‑Jupiter Patched‑Conic Transfer
%======================================================================
% This script performs a first‑order interplanetary “patched‑conic” analysis.
%
% 1. Define Sun, Mercury and Jupiter parameters (μ, orbital radii).
% 2. Compute the heliocentric Hohmann‑type transfer ellipse:
% 3. Determine the circular orbital speeds of the planets (VM, VJ).
% 4. MERCURY DEPARTURE
% 5. JUPITER ARRIVAL
% 6. Total impulsive budget: deltaV = V1 + V2.
% 7. Approximate time‑of‑flight for the Sun‑centred transfer:

% All distances are in km, speeds in km/s, and gravitational parameters in
% km^3/s^2.  The calculations use the classic patched‑conic approximation:
%   Sun‑centred Keplerian arc → planet‑centred hyperbola → Sun‑centred arc.
%======================================================================

us = 1.327*10^11;
RM = 57.909*10^6; 
RJ = 778.479*10^6;  

Vtp = sqrt(2*(-us/(RM+RJ) + us/RM));
Vta = sqrt(2*(-us/(RM+RJ) + us/RJ));

VM = sqrt(us/RM);
VJ = sqrt(us/RJ);

%mercury
ri = 400+2440.5;
um = 0.022032*10^6;
Vi = sqrt(um/ri);

VinfM = Vtp-VM;
Vhp = sqrt(2*(VinfM^2/2 + um/ri));
R_SOI_M = ((0.3301*10^24)/(1988500*10^24))^(2/5) * 57.91*10^6;

V1 = Vhp-Vi;

%jupiter
rf = 10000+71492;
R_SOI_J = (1898.13/1988500)^(2/5) * (778.479*10^6);
uj = 126.687*10^6;
Vf = sqrt(uj/rf);
VinfJ = VJ-Vta;
Vph = sqrt(2*(VinfJ^2/2 + uj/rf));

V2 = Vph-Vf;

deltaV = V1+V2;

%TOF
TOF = 2*pi*sqrt((RM+RJ)^3/us)
%% 2

e = 1.2;
rp = 5380; 
um = 0.042828*10^6;
a = rp/(1-e);
Vinf = sqrt((e-1)*um/rp);
B = acosd(1/e);
turn_angle = asind(1/e)*2


V = sqrt(2*Vinf^2*(1-cosd(turn_angle)))















