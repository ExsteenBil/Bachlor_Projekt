function [Tdiff,Fscreen]=Varmebalance(lufthastighed,Tudenfor,Globalsol,screennummer)
%lufthastighed:Lufthastigheden i hulrummets m/s
%Tudenfor: Temperaturen udenfor i grader Celsius
%Globalsol: Totalstråling på overflade
%screennummer: screennummer

%Tdiff: matematiske temperaturforskel ud fra varmebalance
%Fscreen værdi
%Written by: Ejnar Exsteen & Peter Garlet Rank - S193842 & s193851  2022
%date: 01/06/2022

%Fscreen værdier af screens 
Fscreen = [0.26,0.16,0.13,0.2,0.2,0.18];
Fscreen = Fscreen(screennummer);
%oversigt over screen navne/rækkefølge
screenname = (['816 Original 5%','Silver 4%','Silver 2%','Enviroscreen','878 Original','Clearview']);
%Egenskaber for vindue - Calumen live
alpha1_W = (0.32);
alpha2_W = (0.01);
alpha3_W = (0.01);
%Egenskaber for skærme:
T_s = [0.29,0.06,0.03,0.04,0.03,0.03];
Rho_s = [0.44,0.77,0.82,0.74,0.68,0.74];
alpha_s = 1 - (T_s(screennummer) + Rho_s(screennummer));
epsiv_s = [0.51,0.26 ,0.03 ,0.19,0.27,0.21];


%Vindues emissiviteter:
epsiv11 = 0.84;
epsiv12 = 0.021;
epsiv21 = 0.840;
epsiv22 = 0.840;
epsiv31 = 0.047;
epsiv32 = 0.840;

% emissivitet af bagsiden af skærmen sættes til:
epsiv42 = [0.9, 0.878, 0.878, 0.81, 0.86, 0.85];

%Varmekapacitet:
cp = 1006;

%Karakteristisk længde:
L = 1.5;

%Koefficenter og konstanter:
sigma = 5.6704*10^-8;
g = 9.82;

%Total transmission - Calumen live
tau = (0.24);
lambda_argon = 0.017;

%Geometri af kasse:
w = 1;
h = 1.5;
%Geometri af hulrum:
wg = 0.15;
dp = 0.014;
As = w * h;

%h_ke af Ds418
h_ke = 25;
for i=1:length(lufthastighed)
    %Temperaturer:
Te = 273+Tudenfor(i);
Tsky = 273+Tudenfor(i);
T = 273+Tudenfor(i);
%Hastighed:
v = lufthastighed(i);
%Stråling:
phi_s = Globalsol(i);
phi_si = phi_s * tau;
%Densitet som funktion af temperatur:
rho = (0.029*101325)/(8.314*T);

%Viskositet som funktion af tiden:
nu = (8.65*10^-8)*T-1*10^-5;

%Varmeoverføringsevnen:
lambda = 7*10^-5 * T + 0.0052;

%Prandts tal:
Pr = -0.0001*T + 0.7423;

%Volumen
V = v * w * wg;

%Hydrulisk radius (3-4):
D_h = (2*wg*w)/(wg+w);

%Reynolds som funktion af hastigheden:
Re_v = (v * D_h)/nu;

%Nusselt tal for ventilation:
nusselt = 1/16 * Re_v * (Pr)^(1/3);

%Varmeovergangstal for ventilation:
h_34(i) = ((nusselt*lambda)/D_h);

%Konvektion
syms A B C E F
%kovenktionen
k = [(A+B)/2 (B+C)/2];
%Prandts tal:
Prk = -0.0001 * k + 0.7423;
%konvektion i hulrum
Ra12 = ((g * (1/((A+B)/2)) * (abs(A - B)) * dp^3)/(nu^2)) * Prk(1);

h12 = (((1 + ((0.104 * Ra12^(0.293))/(1 + (6310/Ra12)^1.36))^3 )^(1/3)) * lambda_argon) / (dp);

Ra23 = ((g * (1/((C+B)/2)) * (abs(B - C)) * dp^3)/(nu^2)) * Prk(2);

h23 = (((1 + ((0.104 * Ra23^(0.293))/(1 + (6310/Ra23)^1.36))^3 )^(1/3)) * lambda_argon) / (dp);
%Nussel tal indside :
nusseltki = 91.89730336 * (abs(E-T))^(1/3);
%indre konvektion i kassen
hki = ((nusseltki*lambda)/L);
%Opstilling af ligninger:
eq1 = (epsiv11 * sigma * (A^4 - Tsky^4) * As + h_ke * (A - Te) * As + (1/(1/epsiv12 + 1/epsiv21 - 1)) * sigma * (A^4 - B^4) * As + h12 * (A - B) * As) - (phi_s * alpha1_W * As) == 0;

eq2 = ((1/(1/epsiv12 + 1/epsiv21 - 1)) * sigma * (B^4 - A^4) * As + h12 * (B - A) * As + (1/(1/epsiv22 + 1/epsiv31 - 1)) * sigma * (B^4 - C^4) * As + h23 * (B - C) *As) - (phi_s * alpha2_W * As) == 0;

eq3 = ((1/(1/epsiv31 + 1/epsiv22 - 1)) * sigma * (C^4 - B^4) * As + h23 * (C - B) * As + (1/(1/epsiv32 + 1/epsiv_s(screennummer) - 1)) * sigma * (C^4 - E^4) * As + h_34(i) * ((F - T)/(log((T - C)/(F - C)))) * As) - (phi_s * alpha3_W * As) == 0;

eq4 = ((1/(1/epsiv32 + 1/epsiv_s(screennummer) - 1)) * sigma * (E^4 - C^4) * As + h_34(i) * ((F - T)/(log((T - E)/(F - E)))) * As + epsiv42(screennummer) * sigma * (E^4 - T^4) * As + hki * (E - T) * As) - (phi_si * alpha_s * As) == 0;

eq5 = (h_34(i) * ((F - T)/(log((T - C)/(F - C)))) * As + h_34(i) * ((F - T)/(log((T - E)/(F - E)))) * As) - (V * rho * cp * (F - T)) == 0;

%Vi løser nu de 5 ligninger:

eq = [eq1, eq2, eq3, eq4, eq5];
vars = [A, B, C, E, F];
[A, B, C, E, F] = vpasolve(eq, vars);

Tdiff(i) = F - T;

end
Tdiff = double(Tdiff)';
Tdiff = abs(Tdiff);
end

