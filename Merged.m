clear ; close all; clc

%VARIABLES: 
%alpha_ec = Electric phase angle
%A_lad = area of ladder in metres^2
%A_ps = Active primary slot area
%A_s2 = area of second slot in metres^2
%beta = Damping coefficient
%B_g1 = Air gap flux density
%Bglk = Airgap flux density when s*Ge= 1 (which is the relative slip times the equivalent goodness factor)
%b_sp = Primary slot opening
%b_ss = Secondary slot opening
%b_s1 = Primary slot width
%b_s2 = Secondary slot width
%f_1 = Frequency of the primary
%f_2 = Secondary frequency or secondary slip frequency
%F_xn = Rated thrust
%f_xn = Shear secondary stress
%f_2r = Secondary required frequency
%G_e = Equivalent goodness factor
%G_ei= Effective goodness factor of iron - from equation 17
%g_m = magnetic airgap
%h_sp = Primary tooth height
%h_ss = Secondary tooth height
%h_s1 = Active Primary slot height
%h_s2 = Secondary slot useful height
%I_1 = RMS value of primary phase current / RPS current per phase
%I_1r = Current in Primary winding
%I_1n = RMS phase current for rated thrust (from eq 40)
%j = unknown - should this be j_cor?
%j_cor = Current density in copper winding
%K_w1 = Primary winding factor
    K_w1 = 0.933; %from table at end of thesis
%K_c = Carter coefficient
    K_c = 1+(t_s/(2*pi*g_m))*((1+(b_sp/t_s))*log(1+(b_sp/t_s))+(1-(b_sp/t_s))*log(1-(b_sp/t_s))); %from LIM handbook, should be between 1.1-1.7
    K_c = 1.25; %from table at end of thesis
%K_ss=0.4 Total (Primary and Secondary) core magnetic saturation coefficient
    K_ss = 0.4; %from table at end of thesis
%K_ladder = Ladder coefficient 
    K_ladder = 0.1; %from table at end of thesis
%K_l2  = Secondary leakage inductance coefficient
%K_r = Carter coefficient of rotor
%K_fill = Slot filling factor
    K_fill = 0.6; %according to textbook
%L_m = Magnetization induction
%l_stack = stack width
%l_lad = ladder length, width of ladder, in metres
%l_ec = end-coil length per side
%w_1 = Number of turns per phase in Primary
%L_1l = Primary leakage inductance
%L_2l = Secondary leakage inductance
%lambda_diff1 = Airgap leakage specific permeance in Primary
%lambda_diff2 = Airgap leakage specific permeance in Secondary
%lambda_ecl = Primary end-coil leakage
%lambda_s1 = Slot specific (Nondimensional) permeance in Primary
%lambda_s2 = Slot specific (Nondimensional) permeance in Secondary
%Mu_0 = Permeability of Vacuum or Air
    Mu_0 = 4*pi*10^(-7); %H/m, permeability constant of free space, tube will be at min 0.125 psi though
%N_s2 = Number of slots in secondary per primary length
%omega_1 = primany frequency in radians
%omega_2 = secondary frequency in radians
%p = number of poles
%PF = power factor calculated from EQ42
%Phi_1n = calculated using equation 43
%P_elm = Electromagnetic power from EQ45
%q = q_1 = Number of slots per pole per phase
%R_1 = Primary resistance per phase
%R_2 = secondary resistance
%Rho_copper = resistivity of copper
    Rho_copper = 1.72*10^(-8); %in ohm meters, at 20°C, annealed
%rho_Al = resistivity of aluminium
    rho_Al = 3.125*10.^(-8); %from table at end of thesis
%s = Relative slip
%sigma_Al = density of Aluminium
    sigma_Al = 3.5*10.^7; %from table at end of thesis
%t_s = Slot pitch
    t_s = t/6;
%t_s2 = Secondary slot pitch
%t = Pole pitch
%theta_1m = Rated primary mmf per pole
%u_r = Required speed
%u_s = Synchronous speed from EQ 46
%V_1line = line or supply voltage
%V_1r = Primary voltage
%V_1n = phase voltage
%w_1 = Number of turns per phase in Primary
%w_1r = unknown - should it be w_1?
%y = Coil span

%EQN 1
%Primary Resistance per Phase,
function R_1 = PrimResist(Rho_copper,l_stack,l_ec,j_cor,w_1,I_1r)    
    R_1 = (2*Rho_copper*(l_stack+l_ec)*j_cor*(w_1)^2)/(w_1*I_1r)
end

%EQN 2
%Primary Leakage Inductance,
function L_1l = PrimLeakInduc(Mu_0,p,q_1,lambda_s1,lambda_diff1,l_stack,lambda_ecl,l_ec,w_1)
    L_1l = ((2*Mu_0)/(p*q_1))*((lambda_s1+lambda_diff1)*l_stack+lambda_ecl*l_ec)*w_1
end

%EQN 3
%Primary Slot Specific Permeance (NonDimensional),
function lambda_s1 = PrimSlotSpecPerm(h_s1,beta,h_sp,b_sp,b_s1)
    lambda_s1 = h_s1*(1+3*beta)/(12*b_s1)+(h_sp/b_sp)
end

%EQN 4
%Secondary Slot Specific Permeance (NonDimensional),
function lambda_s2 = SecSlotSpecPerm(h_s2,beta,h_ss,b_ss,b_s2)
    lambda_s2 = h_s2*(1+3*beta)/(12*b_s2)+(h_ss/b_ss)
end

%EQN 5
%Primary Airgap Leakage Specific Permeance,
function lambda_diff1 = PrimAirLeak(K_c,g,b_sp)
    temp1 = K_c*(g/b_sp)
    lambda_diff1 = (5*temp1)/(5+4*temp1)
end

%EQN 6
%Secondary Airgap Leakage Specific Permeance,
function lambda_diff2 = SecAirLeak(K_c,g,b_ss)
    temp2 = K_c*(g/b_ss)
    lambda_diff2 = (5*temp2)/(5+4*temp2)
end

%EQN 7
function lamda_ecl = EndCoilLeakagePermeance()
    lamda_ecl= 0.3*((3*beta)-1)*q1;
end

%EQN 8
function L_m = MagnetizationInductance()
    L_m= h_s2*((6*mu_0*(K_w1*w_1)*t*l_stack)/((pi.^2)*K_c*g*(p*(1+K_ss))));
end

%EQN 9
function r2 = SecondaryResistance()
    R2=12*rho_Al*((K_w1*w_1)/N_s2)*(l_stack/A_s2 + 2*l_lad/A_lad);
end

%EQN 10
function L_2l = SecondaryLeakageInductance()
L_2l=24*Mu_0(l_stack*(lambda_s2+lambda_diff2))*(((K_w1*w_1).^2)/N_s2)*(1+K_ladder);
end

%EQN 12
%electric phase angle
function alpha_ec = ElecPhaseAngle(p,N_s2)
    alpha_ec = 2*pi*p/N_s2
end

%EQN 14
%Secondary slot area
function A_s2 = SecSlotArea(h_s2,b_s2)
    A_s2 = h_s2*b_s2
end

%EQN 15
%Air gap flux density
function B_g1 = FluxDens(Mu_0,theta_1m, g, K_c, K_s, s, G_e)
B_g1 = (Mu_0*theta_1m)/(g*K_c*(1+K_s)*(1+(s^2*G_e^2))^0.5);
end

%EQN 16
function theta_1m = primarymmfpole(K_w1,omega_1,I_1,p)
    theta_1m = (3*sqrt(2)*(K_w1*omega_1)^2*I1)/(pi*p)
end

%EQN 19
function K_l2 = CarterCoefSec2(omega_2,L2l,R2)
    K_l2 = sqrt(1+(omega_2*L2l/R2)^2);
end

%EQN 20
function B_glk = MagFluxDensity()
    B_glk= (Mu_0*theta_1m)/(g*K_c*(1+K_s)*sqrt(1+1.^2));
end

%EQN 21
function K_w1 = PrimaryWindingFactor(y)
    K_w1 = sin((pi/2)*(y/t))*((sin (pi/6))/(q*sin(pi/(6*q))));
end

%EQN 22
function w_1I_1 = CurrentPerTurnPerPhase()
    w_1I_1= (theta_1m*pi*p)/(3*sqrt(2*K_w1)); 
end

%EQN 23
function A_p = PrimaryActiveArea()
    %Ap= Fxn/fxn
    %Ap=2*p*t.^2*(lstack/t)
    A_p = 2*p*t*l_stack ;
end

%EQN 24
function t_s = PrimarySlotPitch()
    t_s=t/6;
end

%EQN 25
function omega_1 = AngFreq()
    omega_2r = 2*pi*f_2r;
    omega_1 = omega_2r;    
end

%EQN 26
function h_s2 = SecondarySlotUsefulHeight(G_ei)
    h_s2 = (G_ei*(pi.^2)*g*K_c*(1+K_s)*K_r*K_l2)/(Mu_0*w_1*(t.^2)*sigma_Al*(1-b_s2/t_s2));
end

%EQN 27
function L_p = PrimLength(t, p)
    L_p = ((2*p)+1)*t;
end

%EQN 28
function A_ps = ActivePrimSlotArea(w_1,I_1,p,q,j_cor,K_fill)
    A_ps = (w_1*I_1)/(p*q*j_cor*K_fill);
end

%EQN 29
function H_s1 = ActivePrimSlotHeight(A_ps,b_s1)
    h_s1 = A_ps/b_s1;
end

%EQN 30
function F_nk = PeakNormalForce(B_glk, Mu_0, p, t, l_stack)
    F_nk = ((B_glk^2)/2*Mu_0)*2*p*t*l_stack;
end

%EQN 31
function F_x = ForwardForce(I_1,R_2,s,G_e,t,f_1)
    F_x = (3*(I_1^2)*R_2*s*G_e)/(2*tau*f_1*(1+(s*G_e)^-2));
end

%EQN 32
function F_xn = ThrustForce(t,w_1,I_1,L_m)
    F_xn = ((3*pi)/(2*t)) * ((w_1*I_1).^2) * (L_m/K_l2) ;
end

%EQN 33
function f_1r = PrimaryRequiredFrequency(f_2,u_r,t)
    f_1r = f_2 + (u_r/(2*t)) ; 
end

%EQN 34
function V_10 = AvailableVoltagePerPhase(V1line)
% This relation describes the RMS value of the voltage per each phase ...
% based on the supply voltage
V_10 = 0.95*(V1line/(sqrt(3))) ; 
end

%EQN35
function b_sp = PrimarySlotOpening(b_s1,g)
    b_sp = (b_s1/g) ;
end

%EQN 38
%Relative slip
function s = RelSlip(f_2, f_1)
s = f_2/f_1;
end

%EQN 39
%Primary voltage
function V_1r = PrimVolt(I_1r, R_1, j, w_1r, L_1l, L_m, w_1, s, R_2, L_2l)
V_1r = I_1r*(R_1 + j*w_1r*L_1l + (j*w_1*L_m* ((R_2/s+j*w_1*L_2l))/(R_2/s+(j*w_1* (L_m+L_2l)))));
end

%EQN 40
function I_1n = CurrentRatedThrust(w_1, I_1)
    I_1n = (w_1*I_1)/w_1;
end

%EQN 41
function P_1n = InputPower(V_1n,I_1n,PF)
    P1n = 3*V_1n*I_1n*PF;
end

%EQN 42
function PF = PowerFactor(Phi_1n)
    PF = cos(Phi_1n);
end

%EQN 43
function Phi_1n = Phi(Z_e)
%Ze is the end effect impedence (calculated through eq 44)
    Phi_1n = arg(Z_e)
end

%EQN 44
function EEI = EndEffectImpedence(R_1,omega_1,L_1l,L_m,R_2,s,L_2l)
    EEI = R_1 + (omega_1*L_1l)*j + ((omega_1*(L_m*((R_2/s)+(omega_1*L_2l)*j)))*j)/((R_2/s)+(omega_1+(L_2l+L_m))*j);
end

%EQN 45
function P_elm = ElectromagPower(I_1n,R_2,s)
    P_elm = 3*(I_1n^2)*R_2*(1-s/s);
end

%EQN 46
function u_s = SynchronousSpeed(t,f_1)
    u_s = 2*t*f_1;
end

%EQN 47
function F_x = MotorThrust(P_elm,u_s)
    F_x = P_elm/u_s;
end
