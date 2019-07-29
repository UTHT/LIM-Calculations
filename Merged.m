clear ; close all; clc
%EQN 1
%Primary Resistance per Phase,
function R_1 = PrimResist(Rho_copper,l_stack,l_ec,j_cor,w_1,I_1r)
    R_1 = (2*Rho_copper*(l_stack+l_ec)*j_cor*(w_1)^2)/(w_1*I_1r)
end
%R1 = Primary resistance per phase.
%Rhocopper = resistivity of copper
%lstack = stack width
%lec = end-coil length per side
%jcor = Current density in copper winding.
%w1 = Number if turns per phase in Primary.
%I1r = Current in Primary winding.

%EQN 2
%Primary Leakage Inductance,
function L_1l = PrimLeakInduc(Mu_0,p,q_1,lambda_s1,lambda_diff1,l_stack,lambda_ecl,l_ec,w_1)
    L_1l = ((2*Mu_0)/(p*q_1))*((lambda_s1+lambda_diff1)*l_stack+lambda_ecl*l_ec)*w_1
end
%L1l = Primary leakage inductance.
%mu0 = Permeability of Vacuum or Air.
%p = number of poles
%lambdas1 = Slot specific (Nondimensional) permeance in Primary.
%lambdadiff1 = Airgap leakage specific permeance in Primary.
%lstack = stack width
%lambdaecl = Primary end-coil leakage.
%lec = Primary end-coil leakage.
%w1 = Number if turns per phase in Primary.

%EQN 3
%Primary Slot Specific Permeance (NonDimensional),
function lambda_s1 = PrimSlotSpecPerm(h_s1,beta,h_sp,b_sp,b_s1)
    lambda_s1 = h_s1*(1+3*beta)/(12*b_s1)+(h_sp/b_sp)
end
%lambdas1 = Slot specific (Nondimensional) permeance in Primary.
%hs1 = Active Primary slot height.
%beta = Damping coefficient.
%hsp = Primary tooth height.
%bsp = Primary slot opening.
%bs1 = Primary slot width.

%EQN 4
%Secondary Slot Specific Permeance (NonDimensional),
function lambda_s2 = SecSlotSpecPerm(h_s2,beta,h_ss,b_ss,b_s2)
    lambda_s2 = h_s2*(1+3*beta)/(12*b_s2)+(h_ss/b_ss)
end
%lambdas2 = Slot specific (Nondimensional) permeance in Secondary
%hs2 = Secondary slot usefull height.
%beta = Damping coefficient.
%hss = Secondary tooth height.
%bss = Secondary slot width.
%bs2 = Secondary slot opening.

%EQN 5
%Primary Airgap Leakage Specific Permeance,
function lambda_diff1 = PrimAirLeak(K_c,g,b_sp)
    temp1 = K_c*(g/b_sp)
    lambda_diff1 = (5*temp1)/(5+4*temp1)
end
%lambdadiff1 = Airgap leakage specific permeance in Primary
%Kc = Carter coefficient for dual slotting.
%g = mechanical airgap
%bsp = Primary slot opening.

%EQN 6
%Secondary Airgap Leakage Specific Permeance,
function lambda_diff2 = SecAirLeak(K_c,g,b_ss)
    temp2 = K_c*(g/b_ss)
    lambda_diff2 = (5*temp2)/(5+4*temp2)
end
%lambdadiff2 = Airgap leakage specific permeance in Secondary
%Kc = Carter coefficient for dual slotting.
%g = mechanical airgap
%bss = Secondary slot opening.

%EQN 7
function lamda_ecl = EndCoilLeakagePermeance()
    %beta is the damping coeffictient wich we would set as well
    q1=  6; % varies depending on our decisions
    beta=0.833;%from equations at end of thesis
    lamda_ecl= 0.3*((3*beta)-1)*q1;
end

%EQN 8
function L_m = MagnetizationInductance()
    h_s2=0.018; %secondary slot height in metres
    mu_0=1.257*10.^(-8) ; %permeability of vacuum or air
    K_w1=0.933 ; %Primary winding factor
    omega_1=617.216 ; %Number if turns per phase in Primary.
    t =0.028 ; %Pole pitch
    l_stack=6.881/1000; %stack width in metres
    g=0.2/1000;% airgap in metres
    p=3;%p= number of poles
    K_c=1.25 ;  %Carter coefficient
    K_ss=0.4  ; %Total (Primary and Secondary) core magnetic saturation coefficient

    L_m= h_s2*((6*mu_0*(K_w1*omega_1)*t*l_stack)/((pi.^2)*K_c*g*(p*(1+K_ss))));
end

%EQN 9
function r2 = SecondaryResistance()
    Kwl= 0.933;%Primary winding factor
    Wl=617.216;%Number if turns per phase in Primary.
    Ns2=40; % Number of slots in secondary per primary length
    Istack = 6.881/1000; %Stack width in metres
    As2= 4.098*10.^(-5); %area of second slot in metres^2
    Ilad=0.004; %ladder length. width of ladder in metres
    Alad= 8.77*10.^(-5);%area of ladder in metres^2
    pAluminium=3.125*10.^(-8); %resistivity of aluminium

    R2=12*pAluminium*((Kwl*Wl)/Ns2)*(Istack/As2 + 2*Ilad/Alad);

%EQN 10
function L_2l = SecondaryLeakageInductance()
u0= 1.257*10.^(-6); %permeability of vacuum or air
Istack= 6.881/1000; % starck width in metres
Lambda_s2= 3.149; %Secondary Slot Specific permeance
Lambda_diff2= 0.15; % Secondary Airgap Leakage specific permeance
K_wl= 0.933; % Primary winding factor
w_l = 617.216; %numbr of turn per phase in primary
N_s2= 40; % Number of slots in secondary per primary length
K_ladder= 0.1; % Ladder coefficient 
L_2l=24*u0(Istack*(Lambda_s2+Lambda_diff2))*(((K_wl*w_l).^2)/N_s2)*(1+K_ladder);

%EQN 12
%electric phase angle
function alpha_ec = ElecPhaseAngle(p,N_s2)
%alphaec = Electric Phase Angle.
%p = Number of poles.
%Ns2 = Number of slots in Secondary per Primary length.
    alpha_ec = 2*pi*p/N_s2
end


%EQN 14
%Secondary slot area
function A_s2 = SecSlotArea(h_s2,b_s2)
%As2 = Area of Secondary slot.
%hs2 = Secondary slot usefull height.
%bs2 =Secondary slot width.
    A_s2 = h_s2*b_s2
end


%EQN 15
%Air gap flux density
function B_g1 = FluxDens(m_0,theta_1m, g, K_c, K_s, s, G_e)
B_g1 = (m_0*theta_1m)/(g*K_c*(1+K_s)*(1+(s^2*G_e^2))^0.5);
end
%B_g1 = Air gap flux density
%m_0 = 1.257*10.^(-6) = permeability of vacuum or air
%theta_1m = rated primary mmf per pole
%g = mechanical air gap
%K_c = 1.25 = Carter coefficient for sual slotting
%K_s = 0.4 = magnetic saturation factor
%s = relative slip
%G_e = equivalent goodness factor


%EQN 16
function theta_1m = primarymmfpole(Kwl,omega1,I1,p)
    theta_1m = (3*sqrt(2)*(Kwl*omega1)^2*I1)/(pi*p)
end

%theta1m = rated primary mmf per pole
%Kwl = primary winding factor
%omegal1 = primany frequency in radians
%I1 = RMS value of primary phase current / RPS current per phase
%p = number of poles


%EQN 19
function K_l2 = CarterCoefSec2(omega2,L2l,R2)
    K_l2 = sqrt(1+(omega2*L2l/R2)^2);
end

%Kl2 = secondary leakage inductance coefficient
%omgea2 = secondary frequency in radians
%L2l = secondary leakage inductance
%R2 = secondary resistance


%EQN 20
function B_glk = MagFluxDensity()
    u_0= 1.257*10.^(-6); %permeability of vacuum or air
    theta_lm =275.722;% Rated primary mmf per pole
    g=0.2/1000;% airgap in metres
    K_c=1.25 ;  %Carter coefficient for sual slotting
    K_s=0.4; % magnetic saturation factor

    B_glk= (u0*thetalm)/(g*Kc*(1+Ks)*sqrt(1+1.^2));
end


%EQN 21
function K_wl = PrimaryWindingFactor(y)
    t=0.028; %Pole pitch
    %y=  coil span
    q= 2; %number of slot per pole per phase
    K_wl = sin((pi/2)*(y/t))*((sin (pi/6))/(q*sin(pi/(6*q))));
end


%EQN 22
function w_1I_1 = CurrentPerTurnPerPhase()
    theta_lm =275.722;% Rated primary mmf per pole
    p=3;%p= number of poles
    K_wl=0.933 ; %Primary winding factor

    w_1I_1= (theta_lm*pi*p)/(3*sqrt(2*K_wl)); 
end


%EQN 23
function A_p = PrimaryActiveArea()
    p=3;%p= number of poles
    t =0.028 ; %Pole pitch
    l_stack=6.881/1000; %stack width in metres
    F_xn=10; %Rated thrust
    f_xn=8.8*10.^3; %Shear secondary stress
    %Ap= Fxn/fxn
    %Ap=2*p*t.^2*(lstack/t)
    A_p=  2*p*t*l_stack ;
end


%EQN 24
function t_s = PrimarySlotPitch()
    t=0.028; %Pole pitch
    t_s=t/6;
end


%EQN 25
function omega_1 = AngFreq()
    f_2r= 4;% Secondary required frequency
    %w1 primary frequency in radians
    omega_2r=2*pi*f_2r;
    omega_1=omega_2r;    
end


%EQN 26
function h_s2 = SecondarySlotUsefulHeight(G_ei)
    %Gei= goodness factor of iron
    g=0.2/1000; % air gap
    K_c= 1.25;% Carter coefficient for dual slotting
    K_s=0.4; % magnetic saturation factor
    K_r=1.5; %Carter coefficient of rotor
    K_l2=1.2; % Secondary leakage inductance coefficient 
    mu_0= 1257*10.^(-6);
    omega_1= % primary frequency in radians
    t=0.028; %Pole pitch
    sigma_Al=3.5*10.^7; %density of Aluminium
    b_s2 = 0.002;% Secondary slot width
    tau_s2= 0.004;%Secondary slot pitch

    h_s2 = (G_ei*(pi.^2)*g*K_c*(1+K_s)*K_r*K_l2)/(mu_0*w_1*(t.^2)*sigma_Al*(1-b_s2/tau_s2));
end


%EQN 27
function L_p = PrimLength(tau, p)
%p is the number of poles
%tau is the pole pitch
    L_p = ((2*p)+1)*tau;
end


%EQN 28
function A_ps = ActivePrimSlotArea(w1,I1,p,q,Jcor,Kfill)
%w1 is the Number of turns per phase in the primary
%I1 is the RMS value of the primary phase current/RPS current per phase
%p is the number of poles
%q is the number of slots per pole per phase
%Jcor is the current density in the copper winding
%Kfill is the slot filling factor (according to the textbook, should be 0.6)
    A_ps = (w_1*I_1)/(p*q*J_cor*K_fill);
end


%EQN 29
function H_s1 = ActivePrimSlotHeight(A_ps,B_s1)
%Aps (from EQ28)
%Bs1 is the primary slot width
    H_s1 = A_ps/B_s1;
end


%EQN 30
function F_nk = PeakNormalForce(B_glk, Mu_0, p, tau, l_stack)
%Bglk is the airgap flux density when s*Ge= 1 (which is the relative slip times the equivalent goodness factor)
%Mu0 is the permeability of the vacuum/air
%p is the number of poles
%tau is the pole pitch
%l is stack width
    F_nk = ((B_glk^2)/2*Mu_0)*2*p*tau*l_stack;
end


%EQN 31
function F_x = ForwardForce(I_1,R_2,s,Ge,tau,f_1)
%I1 is the RMS value of the primary phase current/RPS current per phase
%R2 is secondary resistance
%s is relative slip
%Ge is the equivalent goodness factor
%tau is the pole pitch
%f1 is the frequency of the primary
    F_x = (3*(I_1^2)*R_2*s*Ge)/(2*tau*f_1*(1+(s*Ge)^-2));
end


%EQN 32
function F_xn = ThrustForce(Tau,W_1,I_1,L_m)
    % Tau is the pole pitch
    % W1 is the number of turns per phase in primary
    % I1 is the RMS value of primary phase current
    % Lm is the magnetization induction
    % Kl2 is the secondary leakage inductance coefficient 
    K_l2 = 1.2 ; 
    F_xn = ((3*pi)/(2*Tau)) * ((W_1*I_1).^2) * (L_m/K_l2) ;
end


%EQN 33
function f_1r = PrimaryRequiredFrequency(f_2,U_r,Tau)
    % f2 is the secondary frequency or secondary slip frequency
    % Ur is the required speed
    % Tao is the pole pitch
    f_1r = f_2 + (U_r/(2*Tau)) ; 
end


%EQN 34
function V_10 = AvailableVoltagePerPhase(V1line)

% V1line is the line or supply voltage
% This relation describes the RMS value of the voltage per each phase ...
% based on the supply voltage

V_10 = 0.95*(V1line/(sqrt(3))) ; 

end


%EQN35
function b_sp = PrimarySlotOpening(b_s1,g)
    % bs1 is the primary slot width
    % g is the mechanical airgap

    b_sp = (b_s1/g) ;
end


%EQN 38
%Relative slip
function s = RelSlip(f_2, f_1)
s = f_2/f_1;
end
%f_2 = The secondary frequency or secondary slip frequency
%f_1 = The frequency of the primary


%EQN 39
%Primary voltage
function V_1r = PrimVolt(I_1r, R_1, j, w_1r, L_1l, L_m, w_1, s, R_2, L_2l)
V_1r = I_1r*(R_1 + j*w_1r*L_1l + (j*w_1*L_m* ((R_2/s+j*w_1*L_2l))/(R_2/s+(j*w_1* (L_m+L_2l)))));
end
%V_1r = Primary voltage
%I_1r = Current in primary winding
%R_1 = Primary resistance per phase
%j = Current density
%w_1r = 
%L_1l = Primary leakage inductance
%L_m = Magnetization inductance
%w_1 = Numebr of phases per turn in primary
%s = Relative slip
%R_2 = Secondary resistance
%L_2l = Secondary leakage inductance


%EQN 40
function I_1n = CurrentRatedThrust(w_1, I_1)
%w1 is the Number of turns per phase in the primary
%I1 is the RMS value of the primary phase current/RPS current per phase
    I_1n = (w_1*I_1)/w_1;
end


%EQN 41
function P_1n = InputPower(V_1n,I_1n,PF)
%V1n is the phase voltage
%I1n is the RMS phase current for rated thrust (from eq 40)
%PF is the power factor calculated from EQ42, where the Phi1n is calculated
%from equation 43
    P1n = 3*V_1n*I_1n*PF;
end


%EQN 42
function PF = PowerFactor(Phi_1n)
%Phi1n is calculated using equation 43
    PF = cos(Phi_1n);
end


%EQN 43
function Phi = Phi(Z_e)
%Ze is the end effect impedence (calculated through eq 44)
    Phi = arg(Z_e)
end


%EQN 44
function EEI = EndEffectImpedence(R_1,omega_1,L_1l,L_m,R_2,s,L_2l)
%R1 is the primary resistance per phase
%omega1 is the primary frequency in radians
%L1l is the primary leakage inductance
%Lm is the magnetization inductance
%R2 is the secondary resistance
%s is the relative slip
%L2l is the secondary leakage inductance
    EEI = R_1 + (omega_1*L_1l)j + ((omega_1*(L_m*((R_2/s)+(omega_1*L_2l)j)))j)/((R_2/s)+(omega_1+(L_2l+L_m))j);
end


%EQN 45
function P_elm = ElectromagPower(I_1n,R_2,S)
%I1n is the RMS phase current for rated thrust (from eq 40)
%R2 is the secondary resistance
%S is slip
    P_elm = 3*(I_1n^2)*R_2*(1-S/S);
end


%EQN 46
function SS = SynchronousSpeed(tau,f_1)
%tau is the pole pitch
%f1 is the frequency of the primary
    SS = 2*tau*f_1;
end


%EQN 47
function F_x = MotorThrust(P_elm,Mu)
%Pelm is the electromagnetic power for EQ45
%Mu(s) is the synchronous speed from EQ 46
    F_x = P_elm/Mu;
end

