% two-link arm balanced upright: nonlinear test with torque limit
clear; close all; clc

%% linear model (deviation from upright)
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

%% controllers
p    = [-2+1.95j, -2-1.95j, -6, -7];
Kpp  = place(A, B, p);
Q    = diag([100 100 4 4]);
R    = 0.2*eye(2);
Klqr = lqr(A, B, Q, R);
Kpd  = [45 0 10 0;
         0 10 0 4];
Ks    = {Kpp, Klqr, Kpd};
names = {'Pole placement','LQR','PD'};

%% check nonlinear model against A at upright
disp('nonlinear model at equilibrium:')
disp(arm_nl(zeros(4,1), zeros(2,1))')
e = 1e-6;  Anum = zeros(4);
for i = 1:4
    dx = zeros(4,1);  dx(i) = e;
    Anum(:,i) = (arm_nl(dx,[0;0]) - arm_nl(-dx,[0;0]))/(2*e);
end
fprintf('max |Anum - A| = %.2e\n\n', max(abs(Anum(:)-A(:))));

%% test setup
th0  = 0.3;                          % starting tilt of link 1 (rad)
x0   = [th0; 0; 0; 0];
t    = 0:0.001:6;
band = 0.02*th0;
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
umax = 5;                            % torque limit (N m)
sat  = @(u) max(min(u, umax), -umax);

%% simulate linear (no limit) and nonlinear (limited) for each controller
fprintf('tilt = %.2f rad, torque limit = %.1f N m (nonlinear only)\n', th0, umax);
fprintf('%-16s | %-26s | %-26s\n','', 'linear (Ts, OS, Tpk)', 'nonlinear (Ts, OS, Tpk)');
figure
for k = 1:3
    K = Ks{k};

    xl = initial(ss(A - B*K, zeros(4,2), eye(4), zeros(4,2)), x0, t);
    ul = -(K*xl')';

    [tt, xx] = ode45(@(tt,x) arm_nl(x, sat(-K*x)), [0 t(end)], x0, opts);
    xn = interp1(tt, xx, t);
    un = sat(-(K*xn')');

    [Ts1,OS1,Tp1] = metrics(t, xl(:,1), ul, th0, band);
    [Ts2,OS2,Tp2] = metrics(t, xn(:,1), un, th0, band);
    fprintf('%-16s | %5.2f s %6.2f%% %5.2f N m | %5.2f s %6.2f%% %5.2f N m\n', ...
            names{k}, Ts1,OS1,Tp1, Ts2,OS2,Tp2);

    subplot(3,1,k)
    h = plot(t, xl(:,1), '--', t, xn(:,1)); grid on; hold on
    plot([0 4],[band band],'k:', [0 4],[-band -band],'k:')
    ylabel('\theta_1 (rad)'); title(names{k}); xlim([0 4])
    legend(h, 'linear (no limit)','nonlinear (5 N m limit)')
end
xlabel('time (s)')

%% functions
function xdot = arm_nl(x, u)
    g   = 9.81;
    th1 = pi + x(1);
    th2 = x(2);
    w1  = x(3);  w2 = x(4);
    c2  = cos(th2);  s2 = sin(th2);
    M = [5/3 + c2,     1/3 + 0.5*c2;
         1/3 + 0.5*c2, 1/3];
    h = -0.5*s2;
    C = [ h*(2*w1*w2 + w2^2);
         -h*w1^2];
    G = [g*(1.5*sin(th1) + 0.5*sin(th1+th2));
         0.5*g*sin(th1+th2)];
    acc  = M \ (u - C - G);
    xdot = [w1; w2; acc];
end

function [Ts,OS,Tpk] = metrics(t, th1, u, th0, band)
    idx = find(abs(th1) > band, 1, 'last');
    if isempty(idx)
        Ts = 0;
    elseif idx >= numel(t)
        Ts = NaN;
    else
        Ts = t(idx+1);
    end
    OS  = max(0, -min(th1))/th0*100;
    Tpk = max(abs(u(:)));
end
