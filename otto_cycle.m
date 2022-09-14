clear  
close all
clc

% Inputs
% Engine Parameters
B = 77;                 % Bore(mm)
St = 85.8;              % Stroke(mm)
n_c = 4;                % Number of cylinders()
T_I = 40;               % Inlet air temperature(Deg c)
P_I = 0.918;            % Inlet air pressure(bar)(1*10^5pa)
V_E = 97;               % Volumetric Efficiency(%)
EOS = 6500;             % Engine Operating Speed(rev/min)
r_c = 10.5;             % Compression Ratio()

% Chemical Parameters
P_1 = 101325;                           % Starting Pressure(pa)
T_1 = 298;                              % Starting Temperature(K)
vs = (pi/4)*((B*10^-3)^2)*St*10^-3;     % Swept Volume(m^3)
vc = vs/(r_c-1);                        % Clearance Volume(m^3)
C_p = 1005;                             % Specific heat for constant volume(j/kg)
C_v = 718;                              % Specific heat for constant pressure(j/kg)
R = 287;                                % Gas constant of air(j/kgK)
y = C_p/C_v;                            % Specific heat ratio of air()
afr = 18.5;                             % Air to Fuel ratio()
Q_lhv = 40*10^6;                        % Lower heat value of Fuel(j/kg)

% Calculations
V_1 = vc+vs;
ma = (P_1*V_1)/(T_1*R);
mf = ma/afr;
Q_23 = Q_lhv*mf;
P_1 = 10000;
T_1 = 298;
V_1 = vs+vc;
V_2 = vc;
P_2 = P_1*(V_1/V_2)^y;
T_2 = T_1*(P_2/P_1)^((y-1)/y);
T_3 = (Q_23/(ma*C_v))+T_2;
P_3 = (ma*R*T_3)/vc;
V_3 = vc;
V_4 = vs+vc;
P_4 = P_3*(V_3/V_4)^y;
T_4 = (P_4*V_4)/(R*ma);
P_5 = 1*10^5;
V_5 = vs+vc;
T_5 = (P_5*V_5)/(ma*R);
P_6 = P_5;
V_6 = vc;
T_6 = (P_6*V_6)/(ma*R);
P_7 = P_1;
V_7 = vc;
T_7 = (P_7*V_7)/(ma*R);
T_E = (1-(1/(r_c^(y-1))))*100;      % Thermal Efficiency
S_0 = 6.86305;
S_1 = S_0+(C_p*log((60+273.15)/T_1));
S_2 = S_1;
S_23=S_2+(C_v*log((T_3-T_2)/T_2));          
S_3=S_2+(C_v*log(T_3/T_2));
S_4=S_3;
S_41=S_4+(C_v*log((T_1-T_4)/T_4));


% Plotting
% Matrix Generation
P = [P_1, P_2, P_3, P_4, P_5, P_6, P_7];
V = [V_1, V_2, V_3, V_4, V_5, V_6, V_7];
T = [T_1, T_2, T_3, T_4];
S = [S_1, S_2, S_3, S_4,];
c_r = 1:0.1:20;
Therm_E = (1-(1./(c_r.^(y-1))))*100;

% P-V Diagram
plot(V,P,'-o');
title('P-V, Otto Cycle')
ylabel('Pressure(Bar)')
xlabel('Volume(m^3)')
text(V_1, P_1, '1', 'Fontsize', 12)
text(V_2, P_2, '2', 'Fontsize', 12)
text(V_3, P_3, '3', 'Fontsize', 12)
text(V_4, P_4, '4', 'Fontsize', 12)
text(V_5, P_5, '5', 'Fontsize', 12)
text(V_6, P_6, '6', 'Fontsize', 12)
text(V_7, P_7, '7', 'Fontsize', 12)
Thermal_Efficiency = sprintf('%f',T_E);
hold on 

% Combustion Process
plot([V_2, V_3], [P_2, P_3],'r-')
hold off

% T-S Diagram
figure;
plot(S,T,'-o')
title('T-S, Otto Cycle')
ylabel('Temperature(K)')
xlabel('Entropy(j/K)')
text(S_1, T_1, '1', 'Fontsize', 12)
text(S_2, T_2, '2', 'Fontsize', 12)
text(S_3, T_3, '3', 'Fontsize', 12)
text(S_4, T_4, '4', 'Fontsize', 12)
  
% Thermal Efficiency vs Compression Ratio
figure;
plot(c_r,Therm_E)
% Labelling
title('Thermal efficiency for compression ratio')
ylabel('Thermal Efficiency(%)')
xlabel('Compression Ratio')

