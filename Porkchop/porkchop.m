% ====================================================================
% Earth‑to‑Mars transfer analysis
%======================================================================
% Purpose
%   Generate a pork‑chop (contour) plot that shows how the
%     time‑of‑flight (TOF), departure C₃ and arrival V∞ vary with
%     departure and arrival dates.
%   Accept only physically feasible transfers: 45 ≤ TOF ≤ 500 days.
%   Solve Lambert’s problem for both short‑way (prograde) and
%     long‑way (retrograde) transfers using the universal‑variable
%     routine `lambert`.
%
% What a pork‑chop plot shows
%   A pork‑chop plot is a contour map of interplanetary‑trajectory
%   metrics (TOF, C3, Vinf) plotted over a grid of launch‑date offsets
%   (x‑axis) and arrival‑date offsets (y‑axis).  The contours let a
%   mission designer quickly see:
%     – When the cheapest launch windows occur (low C3).
%     – How launch‑date flexibility trades off against arrival‑speed
%       (V∞) and flight time.
%   The resulting figure is the classic “pork‑chop” used in NASA
%   mission planning.
% --------------------------------------------------------------

mu_sun = 132712440018;
minTOF = 45;
maxTOF = 500;

%Get Julian Date for Departure Date
depart_date = ymdhms2jd(2022,8,1,12,0,0);
depart_date_end = depart_date + maxTOF; %200 day window for departure

%Get Julian Date for Arrival Date;
arrival_date = ymdhms2jd(2023,1,28,12,0,0);
arrival_date_end = arrival_date + maxTOF; %500 day window for arrival

%Days between Aug 1 and Jan28
tspan = 180;

%Initializing Arrays
time = zeros(maxTOF-minTOF+1,maxTOF-minTOF+1); % Adjusted array size
C3_short = zeros(maxTOF-minTOF+1,maxTOF-minTOF+1);
V_inf_short = zeros(maxTOF-minTOF+1,maxTOF-minTOF+1);
C3_long = zeros(maxTOF-minTOF+1,maxTOF-minTOF+1);
V_inf_long = zeros(maxTOF-minTOF+1,maxTOF-minTOF+1);
arrive  = arrival_date;

% Nested Floor Loop
% Between 45 and 500 days for each possible arrival date given the initial
% departure date, then repeat for all possible departure dates.
m = 0;
for depart = depart_date:depart_date_end % Corrected loop range
    m = m + 1;
    n = 0;
    
    for arrive = arrival_date:arrival_date_end
        n = n+1;

        %Initalize TOF
        TOF = arrive - depart; 
        

        %Test TOF for only 45<TOF<500
        if TOF > 45 && TOF <500
            DMs = 1; %if doing short orbit
            DMl = -1; %if doing long orbit
            [Re,Ve] = findEarth(depart);
            [Rm,Vm] = findMars(arrive);

            %Convert TOF to seconds
            TOFs = TOF*86400; 
            
            [Vo,Vf] = LambertSolver(Re,Rm,TOFs,DMs,mu_sun);
            V_inf_e = Vo-Ve;
            C3_short(n,m) = norm(V_inf_e)^2;
            V_inf_short(n,m) = norm(Vm-Vf);
            time(n,m) = TOF;
            
            [Vo,Vf] = LambertSolver(Re,Rm,TOFs,DMl,mu_sun);
            V_inf_e = Vo-Ve;
            C3_long(n,m) = norm(V_inf_e)^2;
            V_inf_long(n,m) = norm(Vm-Vf);
     
        end
    end 
end

%porkchop
figure(1)
hold on
[C0,h0] = contour(time,'k-','ShowText','on');
clabel(C0,h0);
[C1,h1] = contour(C3_short,0:5:50,'r-','ShowText','on');
clabel(C1,h1);
[C2,h2] = contour(V_inf_short,0:10,'b-','ShowText','on');
clabel(C2,h2);
[C3,h3] = contour(C3_long,0:5:50,'r-','ShowText','on');
clabel(C3,h3);
[C4,h4] = contour(V_inf_long,0:10,'b-','ShowText','on');
clabel(C4,h4);
legend('TOF (days)','C3 (km^2/s^2)','V_{\infty} (km/s)','location','northwest');
xlabel('Days past 8/1/2022 (departure offset)');
ylabel('Days past 1/28/2023 (arrival offset)');
title('Contours of TOF, C_3, and V_{\infty} for Earth–Mars transfers (45–500 days)');
