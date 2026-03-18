function jd = ymdhms2jd(year, month, day, hour, minute, second)
% Purpose: Transform a Gregorian date/time (year, month, day, hour, minute, second)
%          into the corresponding Julian Date (JD). Works for all positive JD values.
% Input  : year   – calendar year (e.g., 2022)
%          month  – month number (1‑12)
%          day    – day of month (1‑31)
%          hour   – hour of day (0‑23)
%          minute – minute of hour (0‑59)
%          second – second of minute (0‑59)
% Output : jd – Julian Date (continuous count of days)

day = day + ((second/60 + minute)/60 + hour)/24;

if(month < 3)
  year = year - 1;
  month = month + 12;
end

A = floor(year/100);
B = 2 - A + floor(A/4);

if(year < 1583)
  B = 0;
  if(year == 1582 & month > 9)
    if(day > 14)
      B = 2 - A + floor(A/4);
    end
  end
end
jd = floor(365.25*(year+4716)) + floor(30.6001*(month+1)) + day + B - 1524.5;