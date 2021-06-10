% calculation of poles and zeros
function [out] = poles(A,Alat)

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
Xw = A(1,3) + Q0;
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
Zw = A(3,3);
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

% poles

omega_sp = sqrt(Zw*Mq - Mw*U0);
xi_sp = -(Mq + Zw + Mwdot*U0)/(2*omega_sp);

omega_p = sqrt(-(Zu-Mu*Zw/Mw)*grav/U0);
xi_p = -Xu/(2*omega_p);

syms 's'
pol = det(s*eye(4)-Alat);
coeff = coeffs(pol);
B = coeff(4);
C = coeff(3);
D = coeff(2);

Tr = (B^2+C)/(B+C^2/D); % 1/Tr

Ts = 1/Tr*(-Lr_p + Lb_p/Nb_p*Nr_p)*grav/U0;

omega_dr = sqrt(D/Tr);
xi_dr = (B-Tr)/(2*omega_dr);

% omega_dr = sqrt(Nb_p + Yv*Nr_p);
% xi_dr = (-(Yv + Nr_p) - Lb_p/Nb_p*(Np_p - grav/U0))/(2*omega_dr);
% 
% Tr = -Lp_p + Lb_p/Nb_p*(Np_p-grav/U0); % 1/Tr
% 
% Ts = 1/Tr*(-Lr_p + Lb_p/Nb_p*Nr_p)*grav/U0;

out = double([xi_sp omega_sp; xi_p omega_p; xi_dr omega_dr; 1 abs(Tr); -1 abs(Ts)]);
% if positive real pole, sin(-90°) so xi = -1; if negative real pole,
% sin(90°) so xi = 1;