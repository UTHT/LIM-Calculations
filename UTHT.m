%Initial Specs

V_1line = 9; %Volatge
m_0 = 3; %Number of Phases
F_xn = 10; %Rated thrust
u_n = 1; %Rated Speed
u_r = 1; %Required Speed
l_travel = 1; %Length of travel


%Constants

sigma_Al = 3.2*10^7;  %Density of Al in A^2.s^3/kg.m^3
rho_Copper = 2.3*10^-8;  %Resistivity of Copper in kg.m^3/A^2.s^3 

%Other Specs
g_0 = 0.0002  %Initial airgap (?) in meters
B_glk = 0.7; %Air gap flux density when s_Ge = 1
p = 3;  %Number of poles
tau = ;  %Pole pitch
l_stack = ;  %Stack width
L_StackPerTau = 0.25; %Length of stack per pole pitch
M_0 = 4*pi*10^-7;  %Permeability of vacuum or air in kg.m/A^2.s^2
f_xn = 8.8*10^3;  %Shear secondary strength in Pa
I_1 = ;  %RMs value of primary phase current/RPS value current per phase
R_2 = ;  %Secondary resistance
s = 1; %Relative slip
Ge = 1;  %Equivalent goodness factor
f_1 = ;  %Frequency of primary
f_2r = ;  %Secondary required frequency
j_cor = ;  %Current density in Copper windings


%force Calculations

F_nk = (B_glk^2)*2*p*tau*l_stack/2*M_0   %Peak Normal Force
F_x = 3*(I_1^2)*R_2*s*Ge/2*tau*f_1*(1+(s*Ge^-2))   %Forward Force


%LIM Dims calculations




%Motor Dims

h_s1 = ;  %primary slot height in mm
h_sp = ;  %Primary tooth height in mm
b_s1 = ;  %Primary slot opening in mm
b_t1 = ;  %Primary tooth width in mm
g = ;  %Airgap
h_s2 = ;  %Secondary sheet height
b_s2 = ;  %Secondary sheet width
l_primary = ; %Primary length
l_stak = ;  %Stack width