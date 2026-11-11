% two-link arm balanced upright: pole placement vs LQR vs PD
clear all;
close all; 
clc

%% model (small deviations from upright)
m1 = 1; m2 = 1; L1 = 1; L2 = 1; g = 9.81;
M11 = m1*L1^2/3 + m2*L1^2 + m2*L2^2/3 + m2*L1*L2;
M12 = m2*L2^2/3 + m2*L1*L2/2;
M22 = m2*L2^2/3;
M0  = [M11 M12; M12 M22];
K11 = (m1 + 2*m2)*g*L1/2 + m2*g*L2/2;
K12 = m2*g*L2/2;
K22 = m2*g*L2/2;
Kg  = -[K11 K12; K12 K22];
A = [zeros(2) eye(2); -(M0\Kg) zeros(2)];
B = [zeros(2); inv(M0)];

disp('open-loop eigenvalues:'); disp(eig(A))
fprintf('controllability rank = %d\n\n', rank(ctrb(A,B)));

%% test setup
th0  = 0.1;                  % initial tilt of link 1 (rad)
x0   = [th0; 0; 0; 0];
t    = 0:0.001:6;
band = 0.02*th0;             % 2% settling band

%% pole placement
p   = [-2+1.95j, -2-1.95j, -6, -7];
Kpp = place(A, B, p)
x_pp = initial(ss(A - B*Kpp, zeros(4,2), eye(4), zeros(4,2)), x0, t);
u_pp = -(Kpp*x_pp')';

%% lqr
Q    = diag([100 100 4 4]);
R    = 0.2*eye(2);
Klqr = lqr(A, B, Q, R)
x_lq = initial(ss(A - B*Klqr, zeros(4,2), eye(4), zeros(4,2)), x0, t);
u_lq = -(Klqr*x_lq')';

%% pd (each joint uses its own angle and velocity)
Kp1 = 45; Kp2 = 10; Kd1 = 10; Kd2 = 4;
Kpd = [Kp1 0 Kd1 0;
       0 Kp2 0 Kd2];
disp('PD closed-loop eigenvalues:'); disp(eig(A - B*Kpd))
x_pd = initial(ss(A - B*Kpd, zeros(4,2), eye(4), zeros(4,2)), x0, t);
u_pd = -(Kpd*x_pd')';

%% comparison table
[Ts1,OS1,Tp1] = metrics(t, x_pp(:,1), u_pp, th0, band);
[Ts2,OS2,Tp2] = metrics(t, x_lq(:,1), u_lq, th0, band);
[Ts3,OS3,Tp3] = metrics(t, x_pd(:,1), u_pd, th0, band);
fprintf('\n%-16s %8s %8s %10s\n','Controller','Ts (s)','OS (%)','Peak (N.m)');
fprintf('%-16s %8.2f %8.3f %10.2f\n','Pole placement',Ts1,OS1,Tp1);
fprintf('%-16s %8.2f %8.3f %10.2f\n','LQR',Ts2,OS2,Tp2);
fprintf('%-16s %8.2f %8.3f %10.2f\n\n','PD',Ts3,OS3,Tp3);

%% comparison plots
figure
subplot(3,1,1)
h = plot(t,x_pp(:,1), t,x_lq(:,1), t,x_pd(:,1)); grid on; hold on
plot([0 4],[band band],'k--', [0 4],[-band -band],'k--')
ylabel('\theta_1 (rad)'); xlim([0 4])
legend(h, 'Pole placement','LQR','PD')
subplot(3,1,2)
plot(t,u_pp(:,1), t,u_lq(:,1), t,u_pd(:,1)); grid on
ylabel('\tau_1 (N m)'); xlim([0 4])
subplot(3,1,3)
plot(t,u_pp(:,2), t,u_lq(:,2), t,u_pd(:,2)); grid on
ylabel('\tau_2 (N m)'); xlabel('time (s)'); xlim([0 4])

%% pid demo: constant 0.5 N m disturbance on joint 1
Ki = diag([10 5]);
Cz = [eye(2) zeros(2)];
Apid = [A - B*Kpd, -B*Ki;
        Cz,         zeros(2)];
Bd = [B(:,1); 0; 0];
disp('PID closed-loop eigenvalues:'); disp(eig(Apid))

t2 = 0:0.001:10;
d  = 0.5*ones(size(t2));
xpd  = lsim(ss(A - B*Kpd, B(:,1), eye(4), zeros(4,1)), d, t2, zeros(4,1));
xpid = lsim(ss(Apid, Bd, eye(6), zeros(6,1)), d, t2, zeros(6,1));
fprintf('PD  final theta1 = %.4f rad\n', xpd(end,1));
fprintf('PID final theta1 = %.4f rad\n', xpid(end,1));

figure
plot(t2, xpd(:,1), t2, xpid(:,1)); grid on
xlabel('time (s)'); ylabel('\theta_1 (rad)')
title('constant disturbance: PD vs PID')
legend('PD','PID')

%% settling time, overshoot and peak torque
function [Ts,OS,Tpk] = metrics(t, th1, u, th0, band)
    idx = find(abs(th1) > band, 1, 'last');
    if isempty(idx), Ts = 0; else, Ts = t(min(idx+1, numel(t))); end
    OS  = max(0, -min(th1))/th0*100;
    Tpk = max(abs(u(:)));
end
