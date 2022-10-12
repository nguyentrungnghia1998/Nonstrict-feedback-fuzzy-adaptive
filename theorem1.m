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
theta1 = 10;
theta2 = 105;
%% Initial value
x1{1} = 0.2;
x2{1} = 0;
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
    alpha1{i} = -lamda1*z1{i} - z1{i}/2 - z1{i}*theta1/(2*eta1^2) - k1*ro1;
    z2{i} = x2{i} - alpha1{i};
    if abs(z2{i}) > eps2
        ro2 = sign(z2{i})*abs(z2{i})^r;
    else
        xi2 = 1/2*(3-r)*eps2^(r-1);
        gamma2 = 1/2*(r-1)*eps2^(r-3);
        ro2 = xi2*z2{i} + gamma2*z2{i}^3;
    end
    u{i} = -lamda2*z2{i} - z2{i}/2 - z2{i}*theta2/(2*eta2^2) - k2*ro2;

    if i == length(t)
        break
    end

    % Update
    x1{i+1} = x1{i} + Step*(1 - cos(x1{i}*x2{i}) + (2.5 + 0.5*sin(x1{i}))*x2{i});
    x2{i+1} = x2{i} + Step*(x1{i}^2*exp(x2{i}) + (2+sin(x1{i}*x2{i}))*u{i});
end

x1_m = cell2mat(x1);
x2_m = cell2mat(x2);
figure(1);
plot(t,x1_m,t,x2_m);
z1_m = cell2mat(z1);
z2_m = cell2mat(z2);
figure(2);
plot(t,z1_m,t,z2_m);
yd_m = cell2mat(yd);
figure(3);
plot(t,x1_m,t,yd_m);
figure(4);
u_m = cell2mat(u);
plot(t,u_m);