% calculation of poles and zeros
function  zeros_tf(A,Blon,Blat)

grav = 9.81;
rho_SL = 1.225;
H0 = 3000;
beta_a = 9297;
rho0 = rho_SL*exp(-H0/beta_a);
S = 21.55;
b = 10.48;
c = S/b;
Cmadot = -4.98;

U0 = 158.1484;
% V0 = 0;
W0 = 13.1075;
P0 = -0.0043;
Q0 = 0.0338;
R0 = 0.0397;
Phi0 = 0.7050;
Theta0 = 0.0827;
% Psi0 = 0;

Ixx        = 31481;
Iyy        = 24404;
Izz        = 39361;
Ixz        = 117;

% added terms for new linearization in turn maneuver
Pp = (Q0*Izz*Ixz + Ixz*Q0*(Ixx-Iyy))/(Ixx*Izz-Ixz^2);
Pq = (Izz*(R0*Iyy + Ixz*P0 - Izz*R0) - Ixz*(Ixx*P0 - Ixz*R0 - P0*Iyy))/(Ixx*Izz-Ixz^2);
Pr = (-Q0*Izz^2 + Iyy*Izz*Q0 - Ixz^2*Q0)/(Ixx*Izz-Ixz^2);

Qp = (-R0*Ixx - 2*P0*Ixz + R0*Izz)/Iyy;
Qr = (2*R0*Ixz + P0*Izz - P0*Ixx)/Iyy;

Rp = (Ixx*(Q0*Ixx - Iyy*Q0) + Ixz^2*Q0)/(Ixx*Izz -  Ixz^2);
Rq = (Ixx*(Ixx*P0 - Ixz*R0 - P0*Iyy) + Ixz*(R0*Iyy + Ixz*P0 - Izz*R0))/(Ixx*Izz -  Ixz^2);
Rr = (-Ixx*Ixz*Q0 - Ixz*(Q0*Izz - Iyy*Q0))/(Ixx*Izz -  Ixz^2);

% ------- stability derivatives ------
Xu = A(1,1);  
Xv = A(1,2) - R0;
Xw = A(1,3) + Q0; Xalpha = Xw*U0;
Xp = A(1,4);
Xq = A(1,5) + W0;
Xr = A(1,6);
Xphi = A(1,7);
Xtheta = A(1,8) + grav*cos(Theta0);

Yu= A(2,1)+R0;
Yv= A(2,2); 
Yw= A(2,3)-P0;
Yp= A(2,4)-W0;
Yq= A(2,5);
Yr= A(2,6)+U0;
Yphi= A(2,7)-grav*cos(Theta0)*cos(Phi0);
Ytheta= A(2,8)+grav*sin(Theta0)*sin(Phi0);

Zu = A(3,1) - Q0;
Zv = A(3,2) + P0;
Zw = A(3,3);  Zalpha = Zw*U0;
Zp = A(3,4);
Zq = A(3,5) - U0;
Zr = A(3,6);
Zphi = A(3,7) + grav*cos(Theta0)*sin(Phi0);
Ztheta = A(3,8) + grav*sin(Theta0)*cos(Phi0);

Lu_p= A(4,1); 
Lv_p= A(4,2);   Lb_p = Lv_p*U0;
Lw_p= A(4,3);
Lp_p= A(4,4)-Pp;
Lq_p= A(4,5)-Pq;
Lr_p= A(4,6)-Pr;
Lphi_p= A(4,7);
Ltheta_p= A(4,8);

Mu = A(5,1);
Mv = A(5,2);
Mw = A(5,3);
Mwdot = rho0*S*c^2*Cmadot/4/Iyy;
Mp = A(5,4) - Qp;
Mq = A(5,5);
Mr = A(5,6) - Qr;
Mphi = A(5,7);
Mtheta = A(5,8);

Nu_p = A(6,1);
Nv_p = A(6,2);  Nb_p = Nv_p*U0;
Nw_p = A(6,3);
Np_p = A(6,4) - Rp;
Nq_p = A(6,5) - Rq;
Nr_p = A(6,6) - Rr;
Nphi_p = A(6,7);
Ntheta_p = A(6,8);

% inputs
X_de= Blon(1,1);
X_df= Blon(1,2);
X_dth= Blon(1,3);

Z_de=Blon(2,1);
Z_df=Blon(2,2);
Z_dth=Blon(2,3);

M_de=Blon(3,1);
M_df=Blon(3,2);
M_dth=Blon(3,3);

Y_da = Blat(1,1);
L_da = Blat(2,1);
N_da = Blat(3,1);

Y_dr = Blat(1,2);
L_dr = Blat(2,2);
N_dr = Blat(3,2);

% zeros

% elevator input
Theta_de_sp = -Zw; % not considering Mwdot, very first approx, still good
Tw1_de = -Mq + M_de/Z_de*U0; % good approx of HF zero
Th2_de = sqrt(-Zw*U0*M_de/Z_de + Mw*U0); % good approx (see az_de)
Th3_de = -Th2_de; % good approx
Th1_de = -Xu - (grav-Xalpha)*Zu/Zalpha; % same order of magnitude

% disp('----zeros after elevator input----')
% fprintf(' Theta_de_sp = %f\n Tw1_de = %f\n Th2_de = %f\n Th3_de = %f\n Th1_de = %f\n',Theta_de_sp,Tw1_de,Th2_de,Th3_de,Th1_de)

% thrust input
Tw_dth = -Mq + U0*Mu/Zu; % bad approx
Theta_dth = (Zu*Mw - Mu*Zw)/Mu; % considering Mwdot = 0 % bad approx

% disp('----zeros after thrust input----')
% fprintf(' Tw_dth = %f\n Theta_dth = %f\n',Tw_dth,Theta_dth)

% aileron input
Tbeta2_da = -Lp_p + L_da/N_da*(Np_p - grav/U0); % bad
Tbeta1_da = grav/Tbeta2_da/U0*(L_da/N_da - Lr_p);  %bad

% disp('----zeros after aileron input----')
% fprintf(' Tbeta2_da = %f\n Tbeta1_da = %f\n',Tbeta2_da,Tbeta1_da)

%rudder input
Tbeta3_dr = -(Nr_p + Lp_p + N_dr/Y_dr); % bad approx
Tbeta2_dr = - Lp_p + L_dr/N_dr*(Np_p - grav/U0); % good
Tbeta1_dr = grav/Tbeta2_dr/U0*(L_dr/N_dr*Nr_p - Lr_p); % good

%disp('----zeros after rudder input----')
%fprintf(' Tbeta3_dr = %f\n Tbeta2_dr = %f\n Tbeta1_dr = %f\n',Tbeta3_dr,Tbeta2_dr,Tbeta1_dr)