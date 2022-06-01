%Written by: Ejnar Exsteen & Peter Garlet Rank - S193842 & s193851  2022
%date: 01/06/2022
%main script 
%Out 1-6 figure fra resultater og table af med Fsys og Phivent
clear
clc
format default
close all   
%Bruger inputs
%Tidsintervalg
Tstart = "12:56:00";
Tslut  = "13:00:00";
%Tidsintervalg ja/nej 0=nej,1=ja
n=0;
%Vælg Screen 
screennr=4;
screenname=["816 Original","Silver 4%","Silver 2%","Enviroscreen","878 Original","Clearview"];
%k1=2 fjerne varmebalance 1=brug varmebalance
k1=1;
F1 =("SOL10001.TXT");
F2 =("SOL20001.TXT");
T3 = readtable('Vejrdata.csv','PreserveVariableNames', true);

%script
%se appendix
[TA,hTid,Localtid,Yf,dag]=setu(F1,F2,Tstart,Tslut,n);
%Forsøgs tid i hms
[h,m,s] =hms(Localtid);
%Dtu til Datetime format
Tidmedtu = (T3.("Time(utc)"));
Tidmedtu = datetime(Tidmedtu,'InputFormat','HH:mm:ss d/MM/yyyy');
%intersect Dtu,s vejrdato i forhold til forsøgs dato
dtday = dateshift(Localtid, 'start', 'day');
dt2day = dateshift(Tidmedtu, 'start', 'day');
[dt3, idx] = intersect(dt2day,dtday);
dtday = dateshift(Localtid+days(1), 'start', 'day');
dt2day = dateshift(Tidmedtu, 'start', 'day');
[dt3, idx2] = intersect(dt2day,dtday);
T3 =T3(idx:idx2,:);
%intersect Dtu,s vejrstid i forhold til forsøgs tid i hms
[h2,m2,s2] = hms(Tidmedtu); 
[val1,pos1]=intersect(find(h2==h(1)),find(m2==m(1)));
[val2,pos2]=intersect(find(h2==h(length(h))),find(m2==m(length(m))));
%korrigeret for sommertid
T3 =T3(val1-60:val2-60,:);
%Ude temperaturen
tude = (TA.AAB+TA.ABB)/2;

%Solar Gain
Constant = 9.89*10^(-6);
SolarGain = TA.("Solar Gain")/Constant;

%Script af Toke Rammer Nielsen (DTU)
%Hvori alle parameter er beskrevet i  det vedhæftede script
%se billag
lat =(55.79);
lon = (12.53);
lsm = (15);
Dh = SolarGain*0;
En = Dh;
cloud = 0;
alpha = 0;
beta = 90;
albedo = 0.3;
%Script af Toke Rammer Nielsen (DTU)
[E,ie,alpha_f,gamma_v]=irrad(dag,hTid,Dh,En,cloud,alpha,beta,Yf,lat,lon,lsm,albedo);
%output [alpha_f,gamma_v] er indfaldsvinklen x og y komponent  
dframe =  0;
tframe =  0;
u = 4.7/100;
dx2 = 0;
u2 = 0;
height = 133.2/100;   
width1 = 88/100;
%Script af Toke Rammer Nielsen (DTU)
[F,dz,dx]=nearshadowfactor(alpha_f,gamma_v,dframe,tframe,u,dx2,u2,height,width1);
%F Shadow areal
%Glas areal
Aglas =0.79*1.24;
%G-værdi
Gvalue = 0.26;
%G-60 værdi
Gvalue60 = Gvalue*(1-(tan(((60*pi)/180)/2))^4);
%G-indfalds værdi
Gvaluex = Gvalue*(1-power(tan(((ie*pi)/180)/2),4));
%Dtu diffuse solstråling 
diffusesol =T3.DHI;
%Interpolation af Dtu diffuse solstråling 
diffusesol = interp1(1:length(diffusesol),diffusesol,linspace(1,length(diffusesol),length(SolarGain)))';
%Dtu direkte solstråling 
direktsol = T3.DNI;
%Interpolation af Dtu diffuse solstråling 
SolarDirect = interp1(1:length(direktsol),direktsol,linspace(1,length(direktsol),length(SolarGain)))';
%Direkte solstårling på flade fra Dtu 
SolarDirect2 = (SolarDirect.*cos((ie*pi)/180));
%Difuse stråling
SolarDifuse = (SolarGain-SolarDirect2);
%Tillader Ikke negativ difuse stråling som resultat af difernce mellem
%total stråling på pyrometer og direkte stråling fra Dtu
SolarDifuse(SolarDifuse <= 0) = 0;
%Korrigeret solstråling på flade 
Phiskori = Aglas*(Gvaluex.*SolarDirect2.*F+Gvalue60*SolarDifuse);
%Dtu ude temperatur
tudedtu = T3.air_temperature;
%Interpolation af Dtu ude temperatur
tudedtu = interp1(1:length(tudedtu),tudedtu,linspace(1,length(tudedtu),length(SolarGain)))';
%Antal termoføler 
n = 5;
%Termotårns middel temperatur 
TmeanA = TA.AAB + TA.FA * 25.9 * ((10^3)/(2)/(n));
TmeanB = TA.ABB + TA.FB * 25.9 * ((10^3)/(2)/(n));
%Termperturforskel i termotårne
TinA =  TA.FA*(power(10,3))/n.*(25.9-0.06*TmeanA+2.7*(power(10,-4)).*(power(TmeanA,2))-(power(TmeanA,3))*1*(power(10,-6)));
TinB = TA.FB*(power(10,3))/n.*(25.9-0.06*TmeanB+2.7*(power(10,-4)).*(power(TmeanB,2))-(power(TmeanB,3))*1*(power(10,-6)));
%Tryk
Constant = 50;
Presure = TA.P1;
Presure = Presure*Constant;
%Tryk  til volumenstråm (L/s)
Qreg = -0.0168*(power(Presure,2))+1.5909*Presure+7.0841;
%m^3/s
Qreg = Qreg.*(1/1000);
%korrektion af volumenstrømen iforhold til temperaturen af luften
for e=1:length(tude)
    Rho1(e) = (0.029*101325)/(8.314*(tude(e)+273)); 
    if tude(e)>20
        Qreg(e) = Qreg(e)*(sqrt(1.205/(Rho1(e))));
    else 
        Qreg(e) = Qreg(e)*sqrt((1.293/(Rho1(e))));
    end
end
%areal af tværsnit af hulrum m^2
Acavity = 0.15;
%hastighed i hulrum 
Vcavity = Qreg./Acavity;
%middle hastigheden
Vcavitymean=mean(Vcavity);
%varmebalance script
if k1==1
[Tdiff,Fscreen]=Varmebalance(Vcavity,tude,SolarGain,screennr);
%F-system and % loss heat
%varmekapacitet
Cp = 1006;
%Densitet
Rho = Rho1';
%Varmefjernet af ventilation
PhiRemovedA = Cp .* Rho .* Qreg .* (TinA);
PhiRemovedB = Cp .* Rho .* Qreg .* (TinB);
PhiRemovedMat = (Cp .* Rho .* Qreg .* (Tdiff));
%Varme fjernet fra system med screen
PhimedVentiA = (Phiskori-PhiRemovedA) * Fscreen;
PhimedVentiB = (Phiskori-PhiRemovedB) * Fscreen;
PhimedVentiMat = (Phiskori-PhiRemovedMat) * Fscreen;
%Shading factor
FsysA = (mean(PhimedVentiA./Phiskori))*100;
FsysB = (mean(PhimedVentiB./Phiskori))*100;
FsysMat = (mean(PhimedVentiMat./Phiskori))*100;
FsysAmean = (mean(FsysA))*100;
FsysBmean = (mean(FsysB))*100;
FsysMatmean = (mean(FsysMat))*100;
%Varme fjernet af ventilationen
Phireduktiona = (mean((Phiskori-(Phiskori-PhiRemovedA))./Phiskori))*100;
Phireduktionb = (mean((Phiskori-(Phiskori-PhiRemovedB))./Phiskori))*100;
Phireduktionmat = (mean((Phiskori-(Phiskori-PhiRemovedMat))./Phiskori)')*100;
%Tabel af værdier
Kable = table(FsysA,FsysB,FsysMat,Phireduktiona,Phireduktionb,Phireduktionmat);
disp(Kable)
%Plots
%Temperatur forskel plot
figure 
plot(Localtid,Tdiff,'LineWidth',2,'Color',[1 0 0])
hold on
plot(Localtid,TinA,'LineWidth',2,'Color',[0 1 0])
plot(Localtid,TinB,'LineWidth',2,'Color',[0 0 1])
hold off
ylim([0 4]);
legend(["Diff Mat","Diff A","Diff B"]);
xlabel('Klokkeslæt');
ylabel('Temperaturforskel [K]');
title(strjoin({'Temperaturforskel V =',mat2str(round(Vcavitymean,2)),'[m/s]',' - ',mat2str(screenname(screennr))}));

%Soleffekt og Temperatur forskel plot
figure 
plot(Tdiff,SolarGain,'.','Color',[1 0 0])
hold on
plot(TinA,SolarGain,'.','Color',[0 1 0])
plot(TinB,SolarGain,'.','Color',[0 0 1])
hold off
legend(["Diff Mat","Diff A","Diff B"]);
xlabel('Temperaturforskel [K]');
ylabel('Soleffekt [W/m^2]');
title(strjoin({'Soleffekt og Temperaturforskel V =',mat2str(round(Vcavitymean,2)),'[m/s]',' - ',mat2str(screenname(screennr))}));
else
end
%Gennemsnits temperatur af termo nettet 
g1 = (TA.C1+TA.C4+TA.C7)/3;
g2 = (TA.C2+TA.C5+TA.C8)/3;
g3 = (TA.C3+TA.C6+TA.C9)/3;
%Gennemsnits temperatur i toppen og i bunden
Tbund = (TA.AAB+TA.ABB)/2;
Ttop =  (TA.AAT+TA.ABT)/2;

%Temperatur plot
figure 
plot(Localtid,Ttop,'LineWidth',2,'Color',	[1 0 0])
hold on
plot(Localtid,Tbund,'LineWidth',2,'Color',[0 1 0])
plot(Localtid,g3,'LineWidth',2,'Color',[0 0 1])
plot(Localtid,g2,'LineWidth',2,'Color',[0 1 1])
plot(Localtid,g1,'LineWidth',2,'Color',[1 0 1])
plot(Localtid,tudedtu,'LineWidth',2,'Color',[1 1 0])
ylim([0 26]);
hold off
legend(["Top","Bund","Net Top","Net Midt","Net Bund","Udeluft DTU"]);
xlabel('Klokkeslæt');
ylabel('Absolut temperatur [°C]');
    (strjoin({'Temperatur V =',mat2str(round(Vcavitymean,2)),'[m/s]',' - ',mat2str(screenname(screennr))}));

%Sol effekt plot
figure 
plot(Localtid,SolarGain,'LineWidth',2,'Color',[1 0 0])
hold on
plot(Localtid,SolarDirect2,'LineWidth',2,'Color',[0 1 0])
plot(Localtid,diffusesol,'LineWidth',2,'Color',[0 0 1])
plot(Localtid,SolarDifuse,'LineWidth',2,'Color',[0 1 1])
hold off
ylim([0 1000]);
legend(["Total stråling","Direkte DTU","Diffus DTU","Diffus Forsøg"]);
xlabel('Klokkeslæt');
ylabel('Soleffekt [W/m^2]');
title(strjoin({'Soleffekt på sydvendt lodret flade V =',mat2str(round(Vcavitymean,2)),'[m/s]',' - ',mat2str(screenname(screennr))}));

%Indfaldsvinkel og hastigheds plot
figure
plot(Localtid,ie,'LineWidth',2,'Color',[1 0 0])
hold on
ylabel('Indfaldsvinkel [°]');
yyaxis right
ylabel('Hastighed [m/s]')
plot(Localtid,Vcavity,'LineWidth',2,'Color',[0 1 0])
hold off
legend(["Indfaldsvinkel","Hastighed"]);
xlabel('Klokkeslæt');
title(strjoin({'Indfaldvinkel   &   Hastighed i hulrummet',' - ',mat2str(screenname(screennr))}));

