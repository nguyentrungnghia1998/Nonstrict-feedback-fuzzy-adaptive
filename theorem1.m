%% Time and Step
clc; clear; close all;
Step = 0.005;
T_end = 50;
t = 0:Step:T_end;
data = cell(1,length(t));
%% Variable
alpha1 = data;
u = data;
x1 = data;
x2 = data;
z1 = data;
z2 = data;
yd = data;
theta1 = data;
theta2 = data;
%% Parameter
lamda1 = 5;
eta1 = 1.8;
k1 = 0.6;
r = 0.7;
eps1 = 0.002;
lamda2 = 4;
eta2 = 1.8;
k2 = 2.7;
eps2 = 0.002;
sigma11 = 0.2;
sigma12 = 0.2;
sigma21 = 0.3;
sigma22 = 0.3;
%% Initial value
x1{1} = 0.2;
x2{1} = 0;
theta1{1} = 0.2;
theta2{1} = 0.3;
%% Simulation
for i = 1:length(t)
    yd{i} = sin(t(i));
    z1{i} = x1{i} - yd{i};
    if abs(z1{i}) > eps1
        ro1 = sign(z1{i})*abs(z1{i})^r;
    else
        xi1 = 1/2*(3-r)*eps1^(r-1);
        gamma1 = 1/2*(r-1)*eps1^(r-3);
        ro1 = xi1*z1{i} + gamma1*z1{i}^3;
    end
    alpha1{i} = -lamda1*z1{i} - z1{i}/2 - z1{i}*theta1{i}/(2*eta1^2) - k1*ro1;
    z2{i} = x2{i} - alpha1{i};
    if abs(z2{i}) > eps2
        ro2 = sign(z2{i})*abs(z2{i})^r;
    else
        xi2 = 1/2*(3-r)*eps2^(r-1);
        gamma2 = 1/2*(r-1)*eps2^(r-3);
        ro2 = xi2*z2{i} + gamma2*z2{i}^3;
    end
    u{i} = -lamda2*z2{i} - z2{i}/2 - z2{i}*theta2{i}/(2*eta2^2) - k2*ro2;

    if i == length(t)
        break
    end

    % Update
    x1{i+1} = x1{i} + Step*(1 - cos(x1{i}*x2{i}) + (2.5 + 0.5*sin(x1{i}))*x2{i});
    x2{i+1} = x2{i} + Step*(x1{i}^2*exp(x2{i}) + (2+sin(x1{i}*x2{i}))*u{i});
    theta1{i+1} = theta1{i} + Step*(-sigma11*theta1{i} - sigma12*theta1{i}^r + z1{i}^2/(2*eta1^2));
    theta2{i+1} = theta2{i} + Step*(-sigma21*theta2{i} - sigma22*theta2{i}^r + z2{i}^2/(2*eta2^2));
end

x1_m = cell2mat(x1);
yd_m = cell2mat(yd);
figure(1);
plot(t,x1_m,t,yd_m);
legend('y','y_d');
title('Output response of y');
grid on;
figure(2);
x2_m = cell2mat(x2);
plot(t,x2_m);
legend('x_2');
title('System state x2');
grid on;
theta1_m = cell2mat(theta1);
theta2_m = cell2mat(theta2);
figure(3);
plot(t,theta1_m,'k',t,theta2_m,'g-.');
legend('$\hat{\theta}_{1}$','$\hat{\theta}_{2}$','Interpreter','latex');
title('Adaptive laws $\hat{\theta}_{1}$ and $\hat{\theta}_{2}$','Interpreter','latex');
grid on;
figure(4);
alpha_m = cell2mat(alpha1);
plot(t,alpha_m);
legend('$\alpha_1$','Interpreter','latex');
title('Virtual control law $\alpha_1$','Interpreter','latex');
grid on;
figure(5);
u_m = cell2mat(u);
plot(t,u_m);
legend('u');
title('Control input u');
grid on;