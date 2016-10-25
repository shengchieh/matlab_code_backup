%%
clear all ; close all ; clc;

%%
clear all ; close all ; clc;

N = 50;
r_1 = linspace(0.2, 0.5, N);                        % [0.3,2] 之間取 N 點
theta_1 = linspace(20*pi/180, 80*pi/180, N);       % [20*pi/180,160*pi/180] 之間取 N 點

R = 0.33;               % sensor 間距
M = 0.002;              % dipole 強度
Be = 0.2;               % 地磁強度
theta_e = pi/3;
Bex = Be*cos(theta_e);
Bey = Be*sin(theta_e);
A = Bex / (Be)^2;
B = Bey / (Be)^2;

[X1,Y1] = meshgrid(r_1,theta_1);
X2 = zeros(N,N);   % r_2
Y2 = zeros(N,N);   % theta_2

for i = 1:N
    for j = 1:N
        X2(i,j) = sqrt(X1(i,j)^2 + R^2 - 2*X1(i,j)*R*cos(Y1(i,j)));
        Y2(i,j) = asin((X1(i,j)/X2(i,j))*sin(Y1(i,j)));
    end
end

F = 0*X1;
for i = 1:N
    for j = 1:N
        F(i,j) = ((3*cos(Y1(i,j))*sin(Y1(i,j)))/(X1(i,j)^3) + (3*cos(Y2(i,j))*sin(Y2(i,j)))/(X2(i,j)^3))*A*M + ((3*(cos(Y1(i,j)))^2-1)/(X1(i,j))^3 + (3*(cos(Y2(i,j)))^2-1)/(X2(i,j))^3)*B*M;
    end
end

figure;
surf(X1,Y1,F);
xlabel('r 1') ; ylabel('theta 1')

F2 = 0*F;
hold on;
mesh(X1,Y1,F2);

%%
clear all ; close all ; clc;

r1 = 0.4143;
theta1 = 1.289;
s1 = [r1*sin(theta1) ; r1*cos(theta1)];
R = 0.33;           % sensor 間距
r2 = sqrt( r1^2 + R^2 - 2*r1*R*cos(theta1) );
theta2 = asin( (r1/r2) * sin(theta1) );
s2 = [r2*sin(theta2) ; r2*cos(theta2)];

Be = 0.2;           % 地磁強度
theta_e = pi/3;
Bex = Be*cos(theta_e);
Bey = Be*sin(theta_e);

M = 0.002;              % dipole 強度
dBx1 = M / (r1^3) * 3*cos(theta1)*sin(theta1);
dBy1 = M / (r1^3) * (3*cos(theta1)^2 - 1);
dBx2 = M / (r2^3) * 3*cos(theta2)*sin(theta2);
dBy2 = M / (r2^3) * (3*cos(theta2)^2 - 1);

% total magnetic magnitude
Bx1 = Bex + dBx1;
By1 = Bey + dBy1;
B_1 = [Bx1 ; By1];
Bx2 = Bex + dBx2;
By2 = Bey + dBy2;
B_2 = [Bx2 ; By2];
% trust index
Q = (B_1/Be)'*(B_2/Be) - 1;
% orientation error
phi_e = atan2(Bey,Bex);
% if phi_e < 0
%     phi_e = phi_e + 2*pi;
% end
phi_ed = phi_e*180/pi;

phi_1 = atan2(By1,Bx1);
% if phi_1 < 0
%     phi_1 = phi_1 + 2*pi;
% end
phi_1d = phi_1*180/pi;

phi_2 = atan2(By2,Bx2);
% if phi_2 < 0
%     phi_2 = phi_2 + 2*pi;
% end
phi_2d = phi_2*180/pi;

phi_avg = (phi_1 + phi_2) / 2;
phi_avgd = phi_avg*180/pi;
phi_errd = (phi_avg - phi_e)*180/pi;

% check geometry from figure
figure;
hold on ; axis equal ; grid on;
plot(0,0,'-o','MarkerSize',15,'color','r');
hold on ; axis equal;
% two sensor line
plot(s1(1,1),s1(2,1),'-o','MarkerSize',15,'color','b');
plot(s2(1,1),s2(2,1),'-o','MarkerSize',15,'color','b');
plot([s1(1,1) s2(1,1)],[s1(2,1) s2(2,1)],'LineWidth',1,'Color','b');
% two sensor to magnetic dipole line
plot([0 s1(1,1)],[0 s1(2,1)],'LineWidth',2,'Color','m');
plot([0 s2(1,1)],[0 s2(2,1)],'LineWidth',2,'Color','m');
% delta magnetic dipole vector
plot([s1(1,1) s1(1,1)+dBx1],[s1(2,1) s1(2,1)+dBy1],'LineWidth',2.5,'Color','k');
plot([s2(1,1) s2(1,1)+dBx2],[s2(2,1) s2(2,1)+dBy2],'LineWidth',2.5,'Color','k');
% earth magnetism vector
plot([s1(1,1) s1(1,1)+Bex],[s1(2,1) s1(2,1)+Bey],'LineWidth',2.5,'Color','c');
plot([s2(1,1) s2(1,1)+Bex],[s2(2,1) s2(2,1)+Bey],'LineWidth',2.5,'Color','c');
% total magnetic dipole vector
plot([s1(1,1) s1(1,1)+Bx1],[s1(2,1) s1(2,1)+By1],'LineWidth',2.5,'Color','g');
plot([s2(1,1) s2(1,1)+Bx2],[s2(2,1) s2(2,1)+By2],'LineWidth',2.5,'Color','g');
