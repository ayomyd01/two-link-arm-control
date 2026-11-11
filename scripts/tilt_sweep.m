% largest recoverable tilt with a 5 N m torque limit (nonlinear arm)
clear; close all; clc

%% model and controllers
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

p    = [-2+1.95j, -2-1.95j, -6, -7];
Kpp  = place(A, B, p);
Klqr = lqr(A, B, diag([100 100 4 4]), 0.2*eye(2));
Kpd  = [45 0 10 0; 0 10 0 4];
Ks    = {Kpp, Klqr, Kpd};
names = {'Pole placement','LQR','PD'};

%% sweep settings
umax = 5;
sat  = @(u) max(min(u, umax), -umax);
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
th0s = 0.05:0.01:0.40;               % starting tilts to try (rad)
ok   = false(numel(th0s), 3);

for k = 1:3
    K = Ks{k};
    for i = 1:numel(th0s)
        x0 = [th0s(i); 0; 0; 0];
        [~, xx] = ode45(@(tt,x) arm_nl(x, sat(-K*x)), [0 6], x0, opts);
        ok(i,k) = abs(xx(end,1)) < 0.01 && abs(xx(end,2)) < 0.01;
    end
    nf = find(~ok(:,k), 1);          % first tilt that fails
    if isempty(nf)
        best = th0s(end);
    elseif nf == 1
        best = NaN;
    else
        best = th0s(nf-1);
    end
    fprintf('%-16s largest recovered tilt = %.2f rad (%.1f deg)\n', ...
            names{k}, best, best*180/pi);
end

%% plot
figure
plot(th0s, ok, 'o-'); grid on
xlabel('starting tilt \theta_1 (rad)'); ylabel('recovered = 1, fell = 0')
legend(names); ylim([-0.1 1.1]); title('recovery with a 5 N m torque limit')

%% function
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
