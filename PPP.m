%%
clear all ; close all ; clc;

N = 50;
theta_1 = linspace(0,pi,N);
theta_2 = linspace(0,pi,N);

[X1,Y1] = meshgrid(theta_1,theta_2);

F = 0*X1;
for i = 1:N
    for j = 1:N
        F(i,j) = sin(2*X1(i,j)) / sin(2*Y1(i,j));
        if F(i,j) >= 5
            F(i,j) = 5;
        end
        if F(i,j) <= -5
            F(i,j) = -5;
        end
    end
end

figure;
surf(X1,Y1,F);
xlabel('theta 1 (rad)') ; ylabel('theta 2 (rad)');

F2 = 0*F;
F2 = F2 -1;
hold on;
mesh(X1,Y1,F2);

%%
clear all ; close all ; clc;

N = 100;

theta_1 = linspace(0,pi/2,N);
r = -1                               % r = -(r2/r1)^3
gap = 0.001;
theta_2 = zeros(1,N);
theta_2p = zeros(1,N);

for i=1:N
   theta_2 = 0.5*asin(r*sin(2*theta_1));
   theta_2p = 0.5*asin((r+gap)*sin(2*theta_1));
end

Ap = -trapz(theta_1,theta_2p);
A = -trapz(theta_1,theta_2);
S = (A - Ap)*2;
P = (4*S) / (pi*pi);

plot(theta_1,theta_2);
xlabel('theta 1 (rad)') ; ylabel('theta 2 (rad)');
hold on;
plot(theta_1,theta_2p);

