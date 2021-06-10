clear all
clc
%% Aircraft geometry and aerodynamics;
mass_EO    = 5000;
mass_fuel  = 2500;
mass_res   = 500;
mass0      = mass_EO+mass_fuel+mass_res;
grav       = 9.81;
S          = 21.55;
b          = 10.48;
AR         = b^2/S;
e          = 0.85;
m          = 0.7;
cbar       = 2.05;

Ixx        = 31481;
Iyy        = 24404;
Izz        = 39361;
Ixz        = 117;

Inertia   = [Ixx   0   -Ixz;
             0     Iyy  0;
            -Ixz   0    Izz];

%%% Aerodynamics 

K   = 1/pi/AR/e;
CD0 = 0.023;
Cm0         = 0.0895;
CL0         = 0;
CLa         = 5.16;
CLadot      = 0; 
stab_margin = -0.25;
Cma         = CLa*stab_margin;
CLde        = 0.430;
CLdf        = 1.377;
Cmde        = -1.13;
CLq         = 4.44;
Cmq         = -11.7;
Cmadot      = -4.98;
Cmdf        = 0.4072;

Cy0         = 0;
Cyb         = -0.7678;
Cybdot      = -0.16;
Cyp         = -0.124;
Cyr         = 0.366;
Cyda        = -0.02956;
Cydr        = 0.1158;

Cl0         = 0;
Clb         = -0.06180;
Clbdot      = 0;
Clp         = -0.5045;
Clr         = 0.1695;
Clda        = -0.09917;
Cldr        = 0.006934;

Cn0         = 0;
Cnb         = 0.06719;
Cnbdot      = 0;
Cnp         = -0.1585;
Cnr         = -0.1112;
Cnda        = -0.003872;
Cndr        = -0.08265;

MC_lon = [ CL0 CLa CLadot CLq CLde CLdf;
           Cm0 Cma Cmadot Cmq Cmde Cmdf];
       
MC_lat = [ Cy0 Cyb Cybdot Cyp Cyr Cyda Cydr;
           Cl0 Clb Clbdot Clp Clr Clda Cldr;
           Cn0 Cnb Cnbdot Cnp Cnr Cnda Cndr];

%% Engine

TmaxSL     = 35000;
sfc        = (0.7)/3600;

rho_SL     = 1.225;
beta_a     = 9297;
engine_arm = 0.25;
E_ARM      = [0 0 engine_arm]';
eps0       = (2)*pi/180;

T_TB = [cos(eps0) 0 sin(eps0);
            0     1     0;
       -sin(eps0) 0 cos(eps0)];


%% Flight condition;
% position;
XN0     = 0;
YE0     = 0;
H0      = 3000; 
rho0    = rho_SL*exp(-H0/beta_a);
Tmax0   = TmaxSL*(rho0/rho_SL)^m;
Psi0    = 0;

% Airspeed
CL_Emax = sqrt(CD0/K);
VEmax   = sqrt(2*mass0*grav/S/rho0/CL_Emax);
bV      = 1.31;
V_0      = bV*VEmax;

% Trajectory;
Vz0        = 0;
gamma0     = asin(Vz0/V_0);
Turn_R     = (3)*pi/180;  
Omega_Turn = [0 0 Turn_R]';

%% Trim Procedure - non modificare;

% initial guess;
Phi_guess = atan2(V_0*Turn_R,grav);
n_guess   = 1/cos(Phi_guess);
   
CLtrim    = 2*mass0*grav*n_guess*cos(gamma0)/rho0/S/V_0^2;
NUM_alpha = [CLtrim CLde; 
              -Cm0  Cmde];

NUM_de    = [CLa CLtrim; 
              Cma -Cm0];
          
DEN_trim  = [CLa CLde; 
             Cma Cmde];

alpha_guess = det(NUM_alpha)/det(DEN_trim);
de_guess    = det(NUM_de)/det(DEN_trim);

% iterative
beta_trim = 0;
alpha_i  = alpha_guess;
Phi_i    = Phi_guess;
de_i     = de_guess;
da_i     = 0;
dr_i     = 0;
alpha_i1 = 0;
Phi_i1   = 0;
de_i1    = 0;
Da       = 0;
DPhi     = 0;
Dde      = 0;
J_trim = sqrt((alpha_i1-alpha_i)^2+(de_i1-de_i)^2+(Phi_i1-Phi_i)^2);

while J_trim > 1e-15
alpha_i = alpha_i + 0.5*Da;
Phi_i   = Phi_i + 0.5*DPhi;
de_i    = de_i + 0.5*Dde;
Theta_i = alpha_i+gamma0;

CL_i = CLa*alpha_i + CLde*de_i;
CD_i   = CD0+K*CL_i^2;
L_i = 0.5*rho0*S*V_0^2*CL_i;
D_i = 0.5*rho0*S*V_0^2*CD_i;
n_i = 1/cos(Phi_i);

U_i         = V_0*cos(alpha_i)*cos(beta_trim);
V_i         = V_0*cos(alpha_i)*sin(beta_trim);
W_i         = V_0*sin(alpha_i);
P_i         = -Turn_R*sin(Theta_i);
Q_i         = Turn_R*sin(Phi_i)*cos(Theta_i);
R_i         = Turn_R*cos(Phi_i)*cos(Theta_i);

T_i = (mass0*(Q_i*W_i-R_i*V_i)+mass0*grav*sin(Theta_i)+D_i*cos(alpha_i)-L_i*sin(alpha_i))/cos(eps0);
CmT_i = 2*T_i*engine_arm/rho0/S/V_0^2/cbar;

C_alpha = -CLq*(Q_i*cbar/2/V_0)+2*(-D_i*sin(alpha_i)+mass0*grav*cos(Phi_i)*cos(Theta_i)-T_i*sin(eps0)-mass0*(P_i*V_i-Q_i*U_i))/rho0/S/V_0^2/cos(alpha_i);
C_de    = -Cmq*(Q_i*cbar/2/V_0)-Cm0-CmT_i+2*(P_i*R_i*(Ixx-Ixz)-R_i^2*Ixz+P_i^2*Ixz)/rho0/S/V_0^2;

NUM_alpha_i = [C_alpha CLde; 
              C_de  Cmde];

NUM_de_i    = [CLa C_alpha; 
              Cma C_de];
          
DEN_i  = [CLa CLde; 
          Cma Cmde];

alpha_i1 = det(NUM_alpha_i)/det(DEN_i);
de_i1    = det(NUM_de_i)/det(DEN_i);

Theta_i1 = (alpha_i1+gamma0);


U_i         = V_0*cos(alpha_i1)*cos(beta_trim);
V_i         = V_0*cos(alpha_i1)*sin(beta_trim);
W_i         = V_0*sin(alpha_i1);
P_i         = -Turn_R*sin(Theta_i1);
Q_i         = Turn_R*sin(Phi_i)*cos(Theta_i1);
R_i         = Turn_R*cos(Phi_i)*cos(Theta_i1);

Cl_i       = 2*(Q_i*R_i*(Izz-Iyy)-P_i*Q_i*Ixz)/rho0/S/V_0^2/b;
Cn_i       = 2*(Q_i*P_i*(Iyy-Ixx)+R_i*Q_i*Ixz)/rho0/S/V_0^2/b;

DEN_lat    = [Clda Cldr; 
              Cnda Cndr];
          
NUM_da     = [Cl_i-Clr*(R_i*b/2/V_0)-Clp*(P_i*b/2/V_0) Cldr; 
              Cn_i-Cnr*(R_i*b/2/V_0)-Cnp*(P_i*b/2/V_0) Cndr];
          
NUM_dr    =  [Clda Cl_i-Clr*(R_i*b/2/V_0)-Clp*(P_i*b/2/V_0); 
              Cnda Cn_i-Cnr*(R_i*b/2/V_0)-Cnp*(P_i*b/2/V_0)];

da_i1   = det(NUM_da)/det(DEN_lat);
dr_i1   = det(NUM_dr)/det(DEN_lat);

Cy_i    = Cyp*(P_i*b/2/V_0)+Cyr*(R_i*b/2/V_0)+Cyda*da_i1+Cydr*dr_i1;
Fy_i    = 0.5*V_0^2*S*rho0*Cy_i;
argPhi  = (R_i*U_i-P_i*W_i-Fy_i/mass0)/cos(Theta_i1)/grav;
Phi_i1  = asin(argPhi);

Da   = alpha_i1-alpha_i;
Dde  = de_i1-de_i;
DPhi = Phi_i1-Phi_i;

J_trim = sqrt((alpha_i1-alpha_i)^2+(de_i1-de_i)^2+(Phi_i1-Phi_i)^2);
end

% 
alpha_trim = alpha_i1;
de_trim    = de_i1;
da_trim    = da_i1;
dr_trim    = dr_i1;

Phi0       = Phi_i1;
Theta0     = (alpha_trim+gamma0);

U0         = V_0*cos(alpha_trim)*cos(beta_trim);
V0         = V_0*cos(alpha_trim)*sin(beta_trim);
W0         = V_0*sin(alpha_trim);
P0         = -Turn_R*sin(Theta0);
Q0         = Turn_R*sin(Phi0)*cos(Theta0);
R0         = Turn_R*cos(Phi0)*cos(Theta0);

CL_trim = CLa*alpha_trim+CLq*(Q0*cbar/2/V0)+CLde*de_trim;
L_trim  = 0.5*rho0*S*V_0^2*CL_trim;
D_trim  = 0.5*rho0*S*V_0^2*(CD0+K*CL_trim^2);

T_trim = (mass0*(Q0*W0-R0*V0)+mass0*grav*sin(Theta0)+D_trim*cos(alpha_trim)-L_trim*sin(alpha_trim))/cos(eps0);
dth_trim   = T_i/Tmax0;

%%% variations (only for 'user-performed' linearization);
%%%
du   = 0;
dv   = 0;
dw   = 0;
dp   = 0; %0.01 rad/s
dq   = 0; %0.01 rad/s
dr   = 0; %0.01 rad/s
dde  = 0;
ddf  = 0;
dda  = 0;
ddr  = 0;
ddt  = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

U0         = U0*(1+du);
V0         = V0+dv;
W0         = W0*(1+dw);
P0         = P0+dp;
Q0         = Q0+dq;
R0         = R0+dr;

de_trim    = de_trim*(1+dde);
dth_trim   = dth_trim*(1+ddt);
da_trim    = da_trim+dda;
dr_trim    = dr_trim+ddr;
df0        = 0+ddf;


%Linear Model;
Kq = 0; % used in LINMOD: check poles actually reaching the desired values
Kr = 0; % used in LINMOD: check poles actually reaching the desired values

[A,B,C,D] = linmod('LIN_Complete');
Alon = [A(1,1) A(1,3) A(1,5) A(1,8);
        A(3,1) A(3,3) A(3,5) A(3,8);
        A(5,1) A(5,3) A(5,5) A(5,8);
        A(8,1) A(8,3) A(8,5) A(8,8)];
    
Blon = [B(1,1) B(1,2) B(1,5);
        B(3,1) B(3,2) B(3,5);
        B(5,1) B(5,2) B(5,5);
        B(8,1) B(8,2) B(8,5)];
    
Clon = [1           0       0           0; % u
        0           1       0           0; % w
        0           0       1           0; % q
        0           0       0           1; % theta
        0           1/U0    0           0; % alpha
        0           -1/U0   0           1; % gamma
        A(3,1)-Q0   A(3,3)  A(3,5)-U0   0];% az

Dlon = [zeros(4,3); 0 0 0; 0 0 0; B(3,1) 0 0];
        
        
Alat =  [A(2,2) A(2,4) A(2,6) A(2,7);
         A(4,2) A(4,4) A(4,6) A(4,7);
         A(6,2) A(6,4) A(6,6) A(6,7);
         A(7,2) A(7,4) A(7,6) A(7,7)];
     
Blat = [B(2,3) B(2,4);
        B(4,3) B(4,4);
        B(6,3) B(6,4);
        B(7,3) B(7,4)];
    
Clat = eye(4);

Dlat = zeros(4,2);
         
%% poles and zeros check

[poles] = poles(A(1:8,1:8),Alat);
zeros_tf(A,Blon,Blat)

%% transfer functions

[num_de, den_de]=ss2tf(Alon,Blon,Clon,Dlon,1);
u_de = tf(num_de(1,:),den_de);
w_de = tf(num_de(2,:),den_de);
theta_de = tf(num_de(4,:),den_de);
q_de     = tf(num_de(3,:),den_de);
alpha_de = tf(num_de(5,:),den_de);
gamma_de = tf(num_de(6,:),den_de);
az_de    = tf(num_de(7,:),den_de);

[num_dth, den_dth]=ss2tf(Alon,Blon,Clon,Dlon,3);
u_dth = tf(num_dth(1,:),den_dth);
w_th = tf(num_dth(2,:),den_dth);
theta_dth = tf(num_dth(4,:),den_dth);
q_dth     = tf(num_dth(3,:),den_dth);
alpha_dth = tf(num_dth(5,:),den_dth);
gamma_dth = tf(num_dth(6,:),den_dth);
az_dth    = tf(num_dth(7,:),den_dth);

[num_da, den_da]=ss2tf(Alat,Blat,Clat,Dlat,1);
v_da = tf(num_da(1,:),den_da);
p_da = tf(num_da(2,:),den_da);
r_da = tf(num_da(3,:),den_da);
phi_da = tf(num_da(4,:),den_da);

[num_dr, den_dr]=ss2tf(Alat,Blat,Clat,Dlat,2);
v_dr = tf(num_dr(1,:),den_dr);
p_dr = tf(num_dr(2,:),den_dr);
r_dr = tf(num_dr(3,:),den_dr);
phi_dr = tf(num_dr(4,:),den_dr);

%% ------ SAS ---------

%%%% pitch damper
% Mde = B(5,1);
% xi_sp = 0.24;
% xi_sp_d = 0.5;
% omega_sp = 5.24;
% 
% Kq = 2*omega_sp*(xi_sp_d-xi_sp)/Mde; % improved with sisotool evaluation

%%%% yaw damper
% Np_dr = B(6,4);
% Nb_p = A(6,2)*U0;
% Rr = (-Ixx*Ixz*Q0 - Ixz*(Q0*Izz - Iyy*Q0))/(Ixx*Izz -  Ixz^2);
% Nr_p = A(6,6) - Rr;
% Yv= A(2,2);
% xi_dr_d = 0.2;
% 
% Kr = (Nr_p + Yv + 2*sqrt(Nb_p)*xi_dr_d)/Np_dr;
%

%%%% roll damper
% Tr_d = 1/0.5;
% Pp = (Q0*Izz*Ixz + Ixz*Q0*(Ixx-Iyy))/(Ixx*Izz-Ixz^2);
% Lp_p = A(4,4)-Pp;
% L_da = Blat(2,1);
% 
% Kp = (Tr_d + Lp_p)/L_da;

%% Autopilot Longitudinal

%%%% transfer functions after SAS - Longitudinal Dynamics
Kq = -0.13015; % freezed from sisotool to satisfy handling qualities

Acl1    = Alon-Blon*([0 0 Kq 0 0 0 0; zeros(2,7)]*Clon); % Clon is a 7x4 
% matrix, q is in the 3rd position. The only line which is different from
% zero in the  matrix between Blon and Clon is the first one, i.e. the one
% related to the control of the elevator
Blon1   = Blon;
Clon1   = Clon;
Dlon1   = Dlon;
[num_de1, den_de1] = ss2tf(Acl1,Blon1,Clon1,Dlon1,1);

u_de1       = tf(num_de1(1,:),den_de1);
w_de1       = tf(num_de1(2,:),den_de1);
theta_de1   = tf(num_de1(4,:),den_de1);
q_de1       = tf(num_de1(3,:),den_de1);
alpha_de1   = tf(num_de1(5,:),den_de1);
gamma_de1   = tf(num_de1(6,:),den_de1);
az_de1      = tf(num_de1(7,:),den_de1);

Ktheta = -3.3;
T1 = 8;
T2 = 0.3;
Theta_ref = Theta0;

%% Autopilot Lateral directional

%%%% transfer functions after Kr SAS - Lateral Dynamics
Kr = -0.1655; % freezed from sisotool to satisfy handling qualities

% Acl2    = Alat-Blat*([zeros(1,4); 0 0 Kr 0]*Clat);
% Blat1   = Blat;
% Clat1   = Clat;
% Dlat1   = Dlat;
% [num_dr1, den_dr1] = ss2tf(Acl2,Blat1,Clat1,Dlat1,2);
% 
% v_dr1       = tf(num_dr1(1,:),den_dr1);
% p_dr1       = tf(num_dr1(2,:),den_dr1);
% r_dr1       = tf(num_dr1(3,:),den_dr1); %% just to check of DR damping
% actually increased after SAS
% phi_dr1     = tf(num_dr1(4,:),den_dr1);

%%%% transfer functions after Kp SAS 
Kp = -0.065;
Acl2    = Alat-Blat*([0 Kp 0 0; zeros(1,4)]*Clat);
Blat1   = Blat;
Clat1   = Clat;
Dlat1   = Dlat;
[num_da1, den_da1] = ss2tf(Acl2,Blat1,Clat1,Dlat1,1);

v_da1       = tf(num_da1(1,:),den_da1);
p_da1       = tf(num_da1(2,:),den_da1);
r_da1       = tf(num_da1(3,:),den_da1);
phi_da1     = tf(num_da1(4,:),den_da1);

Kphi = -4.4;
S1 = 0.95;
S2 = 1.2;
Phi_ref = Phi0;
