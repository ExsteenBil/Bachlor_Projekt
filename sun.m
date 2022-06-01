function [LAT,delta,omega,gamma_s,psi_s,alpha_s,tr,ts,omega_s]=sun(day,lst,y,lat,lon,lsm);

%SUN calculates the position of sun on the sky at times corresponding to vectors
%	
%[LAT,delta,omega,gamma_s,psi_s,alpha_s,tr,ts,omega_s]=sun(day,lst,y,lat,lon,lsm)
%	
%LAT: time in LAT (solar time) [h]
%delta: Vector containing declination [deg]
%omega: Vector containing hour angle [deg]
%gamma_s: Vector containing solar elevation [deg]
%psi_s: vector containing solar zenith [deg]
%alpha_s: Vector containing solar azimuth [deg]
%tr: Vector containing sunrise in local standard time [h]
%ts: Vector containing sunset in local standard time [h]
%	 
%day: vector containing day numbers
%lst: vector containing local standard time [h]
%y: year
%lat: lattitude [deg] (optional: default Copemhagen 55.5 deg)
%lon: longitude [deg] (optional: default Copenhagen 12.5 deg)
%lsm: local standard time meridian (optional: default Copenhagen 15 deg)
%	
%Numbers refer to formulaes in The European Radiation Atlas
%
%References: 
%	The European Solar Radiation Atlas
%	École des Mines de Paris, 2000
%	
%   Duffie & Bechmann
%
%Written by: Toke Rammer Nielsen
%Creation date: 6 dec. 2000
%Last modified: 18 sep. 2001

if nargin==3
   lat=55.5;
   lon=12.5;
   lsm=15;
elseif nargin==2
   lat=55.5;
   lon=12.5;
   lsm=15;
   y=inf;
elseif nargin<2
   error('Wrong input')
end

%deg to rad
dg_rd=2*pi/360;
lat=dg_rd*lat; %[rad]
lon=dg_rd*lon; %[rad]
lsm=dg_rd*lsm; %[rad]

%Declination
day_1=2*pi*day/365.25; %(3.3.4)
if y==inf
   delta=23.45*sin(dg_rd*360*(284+day)/365)*dg_rd; %Simple [rad] (Duffie & Bechmann)
else
	%More precise calculation of delta with y as the year
	n0=78.8946+0.2422*(y-1957)-floor((y-1957)/4); %(3.3.2e)
	omega0=2*pi/365.2422; %(3.3.2d)
	t1=-0.5-lon/(2*pi)-n0; %(3.3.2c)
	omegat=omega0*(day+t1); %(3.3.2b)
	delta=0.0064979+0.405906*sin(omegat)+0.0020054*sin(2*omegat)-0.002988*sin(3*omegat)...
   	-0.0132296*cos(omegat)+0.0063809*cos(2*omegat)+0.0003508*cos(3*omegat); %[rad] (3.3.2a)
end

%Time in LAT
ET=-0.128*sin(day_1-0.04887)-0.165*sin(2*day_1+0.34383); %[h] (3.3.22)
c=0; %Correction for summertime not used
LAT=lst+ET+12*((lon-lsm)/pi)-c; %[h] (3.3.23)

%Hour angle
omega=(LAT-12)*pi/12; %[rad] (3.3.1)

%Solar elevation
gamma_s=asin(sin(lat)*sin(delta)+cos(lat)*cos(delta).*cos(omega)); %[rad] (3.3.5)

%Solar zenith
psi_s=pi/2-gamma_s; %[rad] (3.3.6)

%Solar azimuth
sine_alpha_s=cos(delta).*sin(omega)./cos(gamma_s); %(3.3.9)
cosine_alpha_star=(sin(lat)*sin(gamma_s)-sin(delta))./(cos(lat)*cos(gamma_s)); %(3.3.8)
logic=(sine_alpha_s<0);
alpha_s=logic.*(-acos(cosine_alpha_star))+(logic==0).*acos(cosine_alpha_star); %[rad] (3.3.7)

%Sunrise and sunset (with refraction)
check=(sin(-0.0145444)-sin(lat)*sin(delta))./(cos(lat)*cos(delta));
omega_s=acos(check); %(3.3.10)
logic1=check>=1;
logic2=check<=-1;
omega_s=(logic1==0 & logic2==0).*omega_s+logic1*0+logic2*pi; %[rad]

tr=12-12/pi*omega_s-ET-12*((lon-lsm)/pi)+c ; %(3.3.13a) in local standard time
ts=12+12/pi*omega_s-ET-12*((lon-lsm)/pi)+c ; %(3.3.13b) in local standard time

%rad to dg
rd_dg=360/(2*pi);
delta=rd_dg*delta;
omega=rd_dg*omega;
gamma_s=rd_dg*gamma_s;
psi_s=rd_dg*psi_s;
alpha_s=rd_dg*alpha_s;
omega_s=rd_dg*omega_s;