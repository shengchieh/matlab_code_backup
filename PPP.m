%%
clear all ; close all ; clc;

N = 100;
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

figure(1);
surf(X1,Y1,F);
xlabel('theta 1 (rad)') ; ylabel('theta 2 (rad)');

F2 = 0*F;
r = -1;               % r = -(r1/r2)^3
F2 = F2 + r;
hold on;
mesh(X1,Y1,F2);
%
theta_11 = linspace(0,pi/2,N);
gap = 0.001;
theta_22 = zeros(1,N);
theta_22p = zeros(1,N);

for i=1:N
   theta_22 = 0.5*asin((r^-1)*sin(2*theta_11));
   theta_22p = 0.5*asin((r^-1+gap)*sin(2*theta_11));
end

Ap = trapz(theta_11,theta_22p);
A = trapz(theta_11,theta_22);
S = abs(A - Ap)*2;
P = (4*S) / (pi*pi);            % r = -1 ���v���ȵo�� , P = 0.0017

figure(2);
plot(theta_11,theta_22);
axis equal ; grid on;
xlabel('theta 1 (rad)') ; ylabel('theta 2 (rad)');
hold on;
plot(theta_11,theta_22p);

