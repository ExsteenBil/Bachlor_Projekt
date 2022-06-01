function [E,i,alpha_f,gamma_v]=irrad(day,lst,Dh,En,cloud,alpha,beta,y,lat,lon,lsm,albedo);

%IRRAD calculates solar radiation, angle of incidence, horisontal shadow angle
%	and vertical shadow angle on inclined surface at times corresponding to vectors
%	on sloped surfaces
%	
%[E,i,alpha_f,gamma_v]=irrad(day,lst,Dh,En,alpha,beta,lat,lon,lsm,albedo)
%	
%E(1,:): direct solar radiation on inclined surface [W/m2]
%E(2,:): diffuse solar radiation on inclined surface [W/m2]
%E(3,:): reflected solar radiation on inclined surface [W/m2]
%E(4,:): total solar radiation on inclined surface [W/m2]
%i: angle of incidence [deg]
%alpha_f: Horisontal shadow angle [deg]
%gamma_v: Vertical shadow angle [deg]
%
%day: vector containing day numbers
%lst: vector containing local standard time [h]
%Dh: diffuse horisontal irradiation [W/m2]
%En: direct normal irradiation [W/m2]
%cloud: Total cloud cover [00..80] Set to 0 no clouds or 80 if clouded. NOT USED !!
%alpha: azimuth angle of the surface [deg],   
%          east:alpha = -90,   west:alpha = 90
%          south:alpha = 0,    north:alpha = 180
%beta: inclination angle of the surface [deg],
%          horizontal: beta=0, vertical: beta=90
%y: year (set year to inf to use mean calculation of declination)
%lat: lattitude [deg] (optional: default Copenhagen 55.4 deg)
%lon: longitude [deg] (optional: default Copenhagen 12.19 deg)
%lsm: local standard time meridian (optional: default Copenhagen 15 deg)
%albedo: ground reflectivity (optional: default 0.2)
%
%Numbers refer to formulaes in The European Radiation Atlas
%
%References: 
%	The European Solar Radiation Atlas
%	École des Mines de Paris, 2000
%	
%Functions
%	Uses sun.m
%	
%Written by: Toke Rammer Nielsen
%Creation date: 5 dec. 2000
%Last modified: 4 jan. 2001

if nargin==11
   albedo=0.2;
elseif nargin==8
   lat=55.4;
   lon=12.19;
   lsm=15;
   albedo=0.2;
elseif nargin==7
   lat=55.4;
   lon=12.19;
   lsm=15;
   albedo=0.2;
   y=inf;
elseif (nargin<7) | (8<nargin & nargin<11)
   error('Wrong input arguments');
end

logic=cloud>80;
cloud=(logic==0).*cloud+logic.*80;

%Solar position     
[LAT,delta,omega,gamma_s,psi_s,alpha_s,tr,ts,omega_s]=sun(day,lst,y,lat,lon,lsm);

%deg to rad
dg_rd=2*pi/360;
lat=dg_rd*lat; %[rad]
lon=dg_rd*lon; %[rad]
lsm=dg_rd*lsm; %[rad]
alpha=dg_rd*alpha; %[rad]
beta=dg_rd*beta; %[rad]
delta=dg_rd*delta;
omega=dg_rd*omega;
gamma_s=dg_rd*gamma_s;
psi_s=dg_rd*psi_s;
alpha_s=dg_rd*alpha_s;
omega_s=dg_rd*omega_s;

day_1=2*pi*day/365.25; %(3.3.4)

%Incidence angle
cosine_theta_star=cos(omega).*(cos(delta).*(cos(lat)*cos(beta)+cos(alpha)*sin(lat)*sin(beta)))...
   +sin(omega).*(cos(delta)*sin(alpha)*sin(beta))+sin(delta).*(sin(lat)*cos(beta)-cos(lat)*sin(beta)*cos(alpha)); %(3.3.16)
logic=cosine_theta_star>0;
theta=logic.*acos(cosine_theta_star)+(logic==0).*pi/2; %[rad] (3.3.15)
theta2=acos(cosine_theta_star);

%sun path
alpha_f_star=alpha_s-alpha;
logic1=alpha_f_star>pi;
logic2=alpha_f_star<-pi;
alpha_f=(logic1==0 & logic2==0).*alpha_f_star+logic1.*(alpha_f_star-2*pi)+logic2.*(alpha_f_star+2*pi); %[rad] (3.3.20)

gamma_v_star=atan(tan(gamma_s)./cos(alpha_f)); %(3.3.21a)
logic=gamma_v_star<0;
gamma_v=(logic==0).*gamma_v_star+logic.*(gamma_v_star+pi/2); %[rad] (3.3.21b)

%Sky and ground
ri=(1+cos(beta))/2; %(3.3.24)
rg=(1-cos(beta))/2; %(3.3.25)

%%% RADIATION CALCULATIONS %%%
%Direct solar radiation on sloped surface
    logic=(cos(theta)<0 | cos(psi_s)<0);
    E(1,:)=(logic==0).*En.*cos(theta)+logic.*0; %[Wh/m2] (3.6.2)

%Reflected solar radiation on sloped surface
    logic=cos(psi_s)<0;
    Bh=(logic==0).*En.*cos(psi_s)+logic.*0;
    E(3,:)=albedo.*(Bh+Dh)*rg; %[Wh/m2](3.6.4)

%Diffuse irradiation
    if beta==0
        E(2,:)=Dh;
    else
    %Perez "Solar Energy" 1990
    %Air mass
    m=1./(cos(psi_s)+0.15./((93.885-180/pi*psi_s).^1.253));
    logic_m=m>0;
    m=logic_m.*m+(logic_m==0)*0;
        
    %Extraterrestrial radiation
    e_0=1367; %[W/m2]
    E_on=e_0*(1+0.033*cos(dg_rd*360*day/365)); %[W/m2]

	%DETERMINING f11, f12, f13, f21, f22, f23
	logic_Dh=(Dh>0);
	a=1.041*psi_s.^3;
	Dha=Dh+(logic_Dh==0)*1; %To avoid division by zero
	eps=logic_Dh.*(((Dha+En)./Dha+a)./(1+a))+(logic_Dh==0)*1000; %sky clearness

	logic1=(0<=eps & eps<=1.065);
	logic2=(1.065<eps & eps<=1.23);
	logic3=(1.23<eps & eps<=1.5);
	logic4=(1.5<eps & eps<=1.95);
	logic5=(1.95<eps & eps<=2.8);
	logic6=(2.8<eps & eps<=4.5);
	logic7=(4.5<eps & eps<=6.2);
	logic8=(6.2<eps & eps <=inf);
   
    %Perez 1987 "Solar Energy" as in esp-r 
%     logic_Dh=(Dh>0);
% 	Dha=Dh+(logic_Dh==0)*1; %To avoid division by zero
% 	eps=logic_Dh.*(Dha+En)./Dha+(logic_Dh==0)*1000;
%     logic1=(0<=eps & eps<=1.056);
% 	logic2=(1.056<eps & eps<=1.253);
% 	logic3=(1.253<eps & eps<=1.586);
% 	logic4=(1.586<eps & eps<=2.134);
% 	logic5=(2.134<eps & eps<=3.230);
% 	logic6=(3.230<eps & eps<=5.980);
% 	logic7=(5.980<eps & eps<=10.080);
% 	logic8=(10.080<eps & eps <=inf);
%     f11=logic1*(-0.011)+logic2*(-0.038)+logic3*0.166+logic4*0.419+logic5*0.710+logic6*0.857+logic7*0.734...
%     	+logic8*0.421;
% 	f12=logic1*(0.748)+logic2*1.115+logic3*0.909+logic4*0.646+logic5*(0.025)+logic6*(-0.370)+logic7*(-0.073)+...
% 	  logic8*(-0.661);
% 	f13=logic1*(-0.080)+logic2*(-0.109)+logic3*(-0.179)+logic4*(-0.262)+logic5*(-0.290)+logic6*(-0.279)+logic7*(-0.228)+...
%  	  logic8*0.097;
%     f21=logic1*(-0.048)+logic2*(-0.023)+logic3*0.062+logic4*0.140+logic5*0.243+logic6*0.267+logic7*0.231...
%       +logic8*0.119;
% 	f22=logic1*(0.073)+logic2*0.106+logic3*(-0.021)+logic4*(-0.167)+logic5*(-0.511)+logic6*(-0.792)+logic7*(-1.180)+...
%     	logic8*(-2.125);
% 	f23=logic1*(-0.024)+logic2*(-0.037)+logic3*(-0.050)+logic4*(-0.042)+logic5*(-0.004)+logic6*(0.076)+logic7*0.199+...
%     	logic8*0.446;
%     hp=0.5*pi;
%     halfangle=25*pi/180;
%     hpmalf=hp-halfangle;
%     hppalf=hp+halfangle;
%     psic=((hppalf-theta)/halfangle)/2;
%     logic_zet=psi_s>hpmalf;
%     %psih=logic_zet.*(hppalf-psi_s)/(2*halfangle)+(logic_zet==0).*1; %Perez 1987
%     psih=logic_zet.*2.*(hppalf-psi_s)/(halfangle)+(logic_zet==0).*1; %esp-r formulation
%     logic_value1=(theta>=hpmalf).*(theta<=hppalf);
%     logic_value2=(theta<hpmalf);
%     xic=logic_value1.*psih.*psic.*sin(psic.*halfangle)+logic_value2.*psih.*cos(theta)+(logic_value1==0).*(logic_value2==0)*0;
%     logic_value3=(psi_s<hpmalf);
%     xih=logic_value3.*cos(psi_s)+(logic_value3==0).*psih.*sin(psih.*halfangle);
%     a=xic;
%     b=xih;
    
    %Perez fra solstråling
	%f11=logic1*(-0.196)+logic2*0.236+logic3*0.454+logic4*0.866+logic5*1.026+logic6*0.978+logic7*0.748...
    % 	+logic8*0.318;
	%f12=logic1*(1.084)+logic2*0.519+logic3*0.321+logic4*(-0.381)+logic5*(-0.711)+logic6*(-0.986)+logic7*(-0.913)+...
	%   logic8*(-0.757);
	%f13=logic1*(-0.006)+logic2*(-0.18)+logic3*(-0.255)+logic4*(-0.375)+logic5*(-0.426)+logic6*(-0.35)+logic7*(-0.236)+...
    %	logic8*0.103;
    %f21=logic1*(-0.114)+logic2*(-0.011)+logic3*0.072+logic4*0.203+logic5*0.273+logic6*0.28+logic7*0.173...
    %   +logic8*0.062;
	%f22=logic1*(0.18)+logic2*0.02+logic3*(-0.098)+logic4*(-0.403)+logic5*(-0.602)+logic6*(-0.915)+logic7*(-1.045)+...
    %	logic8*(-1.698);
	%f23=logic1*(-0.019)+logic2*(-0.038)+logic3*(-0.046)+logic4*(-0.049)+logic5*(-0.061)+logic6*(-0.024)+logic7*0.065+...
    %	logic8*0.236;
   
    %Perez fra 1990 "Solar Energy"
    f11=logic1*(-0.008)+logic2*0.130+logic3*0.330+logic4*0.568+logic5*0.873+logic6*1.132+logic7*1.060...
    	+logic8*0.678;
	f12=logic1*(0.588)+logic2*0.683+logic3*0.487+logic4*(0.187)+logic5*(-0.392)+logic6*(-1.237)+logic7*(-1.6)+...
	    logic8*(-0.327);
	f13=logic1*(-0.062)+logic2*(-0.151)+logic3*(-0.221)+logic4*(-0.295)+logic5*(-0.362)+logic6*(-0.412)+logic7*(-0.359)+...
  	    logic8*(-0.25);
    f21=logic1*(-0.06)+logic2*(-0.019)+logic3*0.055+logic4*0.109+logic5*0.226+logic6*0.288+logic7*0.264...
        +logic8*0.156;
	f22=logic1*(0.072)+logic2*0.066+logic3*(-0.064)+logic4*(-0.152)+logic5*(-0.462)+logic6*(-0.823)+logic7*(-1.127)+...
    	logic8*(-1.377);
	f23=logic1*(-0.022)+logic2*(-0.029)+logic3*(-0.026)+logic4*(-0.014)+logic5*(0.001)+logic6*(0.056)+logic7*0.131+...
   	    logic8*0.251;
       
    %Determining F1, F2
	%logic_Delta=(0<=psi_s & psi_s<=pi/2);
	%Delta=logic_Delta.*m.*Dh./E_on+(logic_Delta==0).*0;
    Delta=m.*Dh./E_on;
    
	F1=f11+f12.*Delta+psi_s.*f13;
	F2=f21+f22.*Delta+psi_s.*f23;

	a=max(0,cos(theta));
	b=max(0.087,cos(psi_s));

	%Diffuse solar radiation on inclined surface
	E(2,:)=Dh.*((1-F1).*(1+cos(beta))/2+F1.*a./b+F2.*sin(beta)); %[W/m2]
    E(2,:)=(E(2,:)>0).*E(2,:);
    end
%Global radiation on sloped surface
    E(4,:)=E(1,:)+E(2,:)+E(3,:); %[Wh/m2] (3.6.1)

%Angles in degrees
i=180/pi*theta; %[deg]
alpha_f=180/pi*alpha_f; %[deg]
gamma_v=180/pi*gamma_v; %[deg]







