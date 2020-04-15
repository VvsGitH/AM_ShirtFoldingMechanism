%% TEST PARAMETERS

% PLATE
L = 0.35;
H = 0.04;
M = 4.5;
rp = 0.5*sqrt(L^2+H^2);
Ip = M*(L^2 + H^2)/3;

% SLEEVE
l = 0.30;
rm = 0.5*l^2;

% MACHANISM
g = 9.81;
mu = 0.35;
bL = L/2;
bH = H/2;
Cm0 = M*g*bL;

%% GENERATION OF THE POLYNOMIAL CURVE

% CONDITIONS
t0 = 0;
t1 = 0.4;
tc = t1/2;
theta0 = 0;
theta1 = pi;
omegac = sqrt(g/rm); % from adherence equations

% COEFFICIENTS
syms a0 a1 a2 a3 a4 a5;

% SYSTEM OF CONSTRAINTS
eq1 = a0 + a1*t0 + a2*t0^2 + a3*t0^3 + a4*t0^4 + a5*t0^5 == theta0;
eq2 = a0 + a1*t1 + a2*t1^2 + a3*t1^3 + a4*t1^4 + a5*t1^5 == theta1;
eq3 = 0*a0 + a1 + 2*a2*t0 + 3*a3*t0^2 + 4*a4*t0^3 + 5*a5*t0^4 == 0;
eq4 = 0*a0 + a1 + 2*a2*tc + 3*a3*tc^2 + 4*a4*tc^3 + 5*a5*tc^4 == omegac;
eq5 = 0*a0 + 0*a1 + 2*a2 + 6*a3*t1 + 12*a4*t1^2 + 20*a5*t1^3 == 0;
eq6 = 0*a0 + 0*a1 + 2*a2 + 6*a3*tc + 12*a4*tc^2 + 20*a5*tc^3 == 0;

% SYSTEM SOLVING
[A, B] = equationsToMatrix([eq1,eq2,eq3,eq4,eq5,eq6],[a0, a1, a2, a3, a4, a5]);
X = linsolve(A, B);
x = double(X);

%% CURVE FUNCTIONS

t = (0:0.01:t1)';
theta = x(1) + x(2)*t + x(3)*t.^2 + x(4)*t.^3 + x(5)*t.^4 + x(6)*t.^5;
omega = x(2) + 2*x(3)*t + 3*x(4)*t.^2 + 4*x(5)*t.^3 + 5*x(6)*t.^4;
alpha = 2*x(3) + 6*x(4)*t + 12*x(5)*t.^2 + 20*x(6)*t.^3;

% ADHERENCE CONDITIONS
omega_min = real(sqrt(g/rm*sin(theta) - mu*g/rm*cos(theta) - mu*alpha));
omega_max = real(sqrt(g/rm*sin(theta) + mu*g/rm*cos(theta) + mu*alpha));

%% TORQUE CALCULATION

Cm = Ip*alpha + bL*M*g*cos(theta) + bH*M*g*sin(theta);
Cm_max = max(Cm);
Cm_min = min(Cm);
str = ['Coppia Massima = ',num2str(Cm_max,5), ' Nm'];
disp(str)
str = ['Coppia Minima = ',num2str(Cm_min,5), ' Nm'];
disp(str)

%% PLOTTING

% THETA OMEGA TORQUE
figure(1)
subplot(3,1,1)
    plot(t,theta*180/pi)
    grid on
    xlabel('TIME [s]'), ylabel('THETA [grad]')
    legend('THETA','Location','southEast');
subplot(3,1,2)
    plot(t,omega*180/pi,t,omega_min*180/pi,'g',t,omega_max*180/pi,'r')
    grid on
    xlabel('TIME [s]'), ylabel('OMEGA [grad/s]')
    legend('OMEGA', 'OMEGA_M_I_N [adh]', 'OMEGA_M_A_X [adh]','Location','south')
subplot(3,1,3)
    plot(t,Cm)
    grid on
    xlabel('TIME [s]'), ylabel('Cm [Nm]')
    legend('Cm','Location','northEast');

% THETA TO SAMPLE
figure(2)
    plot(t,theta*180/pi)
    grid on
    xlabel('TIME [s]'), ylabel('THETA [grad]')
    legend('THETA','Location','southEast');