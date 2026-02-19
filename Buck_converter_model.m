clear all;
close all;
clc

tsim = 0.7;
%% Parameters of the Buck Converter that we shall use and simulate the environment 

%Without Parasitics

Vin = 50;
L = 100e-6;          %Inductor in Henry 
C = 100e-6;        %Capacitor in Farad
R = 10;           %Resistance in ohm
i_L_int =0;        %Initial condition of inductor current 
v_c_int = 0;       %Initial condition of capacitor voltage

%With Parasitics 


L_p = L;           %Inductor in Henry 
rl = 10e-3;        %Inductor parasitics
r1 = 5e-3;         %FET parasitic resistance
rd = 4e-3;         %Diode Parasitic resistance
vd = 0.7;          %Diode Parasitic voltage
C_p = C;           %Capacitor in Farad
rc = 5e-3;         %ESR

i_L_int =0;        %Initial condition of inductor current 
v_c_int = 0;       %Initial condition of capacitor voltage

%% Control Parameters of Buck Converter 

fsw = 100e3;
D = 0.5;



%% Open the Simulink model of the Buck Converter 
open("Simulink_buck_model.slx");