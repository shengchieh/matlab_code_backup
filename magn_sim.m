%%
clear all ; close all ; clc;

%%
clear all ; close all ; clc;

% given from random determine
rm = 0.5;                       % magnetic dipole center to sensor distance
theta_m = pi/6;                 % magnetic dipole center to sensor distance
theta_d = pi/3;                 % magnetic dipole orientation
theta_e = pi/3;                 % earth magnetism orientation
% start random determine
ran_o = 2*pi*rand();
rxm = 0.25 + 0.3*rand();
rym = 0.25 + 0.5*rand();
rm = sqrt( rxm^2 + rym^2 );
theta_m = ran_o;
theta_d = ran_o;
theta_e = ran_o;
% given from measurement
mag = 0.01;                     % magnetic dipole
Be = 0.25;                      % earth magnetism
% sensor coordinates
s1 = [0 ; 0];
s2 = [-0.11 ; -0.34];
% magnetic dipole coordinates
xm = rm*cos(theta_m);
ym = rm*sin(theta_m);
m = [xm ; ym];
% magnetic dipole orientation
chk = 0.07;
m2 = [xm + chk*sin(theta_d) ; ym + chk*cos(theta_d)];
% distance from magnetic dipole to sensor
r1 = sqrt( ( m(1,1)-s1(1,1) )^2 + ( m(2,1)-s1(2,1) )^2 );
r2 = sqrt( ( m(1,1)-s2(1,1) )^2 + ( m(2,1)-s2(2,1) )^2 );
% orientation from magnetic dipole to sensor
theta_11 = atan2( s1(2,1)-m(2,1) , s1(1,1)-m(1,1) );
if theta_11 < 0
    theta_11 = theta_11 + 2*pi;
end
theta_22 = atan2( s2(2,1)-m(2,1) , s2(1,1)-m(1,1) );
if theta_22 < 0
    theta_22 = theta_22 + 2*pi;
end
theta_1 = (2*pi - theta_11) + (pi/2 - theta_d);
theta_2 = (2*pi - theta_22) + (pi/2 - theta_d);
theta_1d = theta_1*180/pi;
theta_2d = theta_2*180/pi;
% magnetic dipole magnitude : Let mag = u0*M/4*pi
dBx1 = mag / (r1^3) * 3*cos(theta_1)*sin(theta_1);
dBy1 = mag / (r1^3) * (3*cos(theta_1)^2 - 1);
dBx2 = mag / (r2^3) * 3*cos(theta_2)*sin(theta_2);
dBy2 = mag / (r2^3) * (3*cos(theta_2)^2 - 1);
% earth magnetism
Bex = Be*cos(theta_e);
Bey = Be*sin(theta_e);
% total magnetic magnitude
Bx1 = Bex + dBx1;
By1 = Bey + dBy1;
B_1 = [Bx1 ; By1];
Bx2 = Bex + dBx2;
By2 = Bey + dBy2;
B_2 = [Bx2 ; By2];
% trust index
Q = (B_1/Be)'*(B_2/Be);
V = Q - 1;
% orientation error
phi_e = atan2(Bey,Bex);
phi_1 = atan2(By1,Bx1);
phi_2 = atan2(By2,Bx2);
phi_avg = (phi_1 + phi_2) / 2;
phi_errd = (phi_avg - phi_e)*180/pi;
% check geometry from figure
figure;
plot(s1(1,1),s1(2,1),'-o','MarkerSize',15,'color','b');
hold on ; grid on ; axis equal;
% two sensor line
plot(s2(1,1),s2(2,1),'-o','MarkerSize',15,'color','b');
plot([s1(1,1) s2(1,1)],[s1(2,1) s2(2,1)],'LineWidth',1,'Color','b');
% magnetic dipole direction
plot(m(1,1),m(2,1),'-o','MarkerSize',15,'color','r');
plot(m2(1,1),m2(2,1),'-o','MarkerSize',10,'color','r');
plot([m(1,1) m2(1,1)],[m(2,1) m2(2,1)],'LineWidth',1,'Color','r');
% two sensor to magnetic dipole line
plot([m(1,1) s1(1,1)],[m(2,1) s1(2,1)],'LineWidth',1.3,'Color','m');
plot([m(1,1) s2(1,1)],[m(2,1) s2(2,1)],'LineWidth',1.3,'Color','m');
% delta magnetic dipole vector
plot([s1(1,1) s1(1,1)+dBx1],[s1(2,1) s1(2,1)+dBy1],'-.','LineWidth',1.5,'Color','k');
plot([s2(1,1) s2(1,1)+dBx2],[s2(2,1) s2(2,1)+dBy2],'-.','LineWidth',1.5,'Color','k');
% earth magnetism vector
plot([s1(1,1) s1(1,1)+Bex],[s1(2,1) s1(2,1)+Bey],'-.','LineWidth',1.5,'Color','c');
plot([s2(1,1) s2(1,1)+Bex],[s2(2,1) s2(2,1)+Bey],'-.','LineWidth',1.5,'Color','c');
% total magnetic dipole vector
plot([s1(1,1) s1(1,1)+Bx1],[s1(2,1) s1(2,1)+By1],'LineWidth',1.5,'Color','g');
plot([s2(1,1) s2(1,1)+Bx2],[s2(2,1) s2(2,1)+By2],'LineWidth',1.5,'Color','g');

%%
clear all ; close all ; clc;

% given from random determine
rm = 0.5;                       % magnetic dipole center to sensor distance
theta_m = pi/6;                 % magnetic dipole center to sensor distance
theta_d = pi/3;                 % magnetic dipole orientation
theta_e = pi/3;                 % earth magnetism orientation
% start random determine
N = 1000;
R = zeros(N,6);
for i = 1 : N
ran_o = 2*pi*rand();
rxm = 0.45 + 0.3*rand();
rym = 0.45 + 0.5*rand();
rm = sqrt( rxm^2 + rym^2 );
theta_m = ran_o;
theta_d = ran_o;
theta_e = ran_o;
% given from measurement
mag = 0.01;                     % magnetic dipole
Be = 0.25;                      % earth magnetism
% sensor coordinates
s1 = [0 ; 0];
s2 = [-0.11 ; -0.34];
% magnetic dipole coordinates
xm = rm*cos(theta_m);
ym = rm*sin(theta_m);
m = [xm ; ym];
% magnetic dipole orientation
chk = 0.07;
m2 = [xm + chk*sin(theta_d) ; ym + chk*cos(theta_d)];
% distance from magnetic dipole to sensor
r1 = sqrt( ( m(1,1)-s1(1,1) )^2 + ( m(2,1)-s1(2,1) )^2 );
r2 = sqrt( ( m(1,1)-s2(1,1) )^2 + ( m(2,1)-s2(2,1) )^2 );
% orientation from magnetic dipole to sensor
theta_11 = atan2( s1(2,1)-m(2,1) , s1(1,1)-m(1,1) );
if theta_11 < 0
    theta_11 = theta_11 + 2*pi;
end
theta_22 = atan2( s2(2,1)-m(2,1) , s2(1,1)-m(1,1) );
if theta_22 < 0
    theta_22 = theta_22 + 2*pi;
end
theta_1 = (2*pi - theta_11) + (pi/2 - theta_d);
theta_2 = (2*pi - theta_22) + (pi/2 - theta_d);
theta_1d = theta_1*180/pi;
theta_2d = theta_2*180/pi;
% magnetic dipole magnitude : Let mag = u0*M/4*pi
dBx1 = mag / (r1^3) * 3*cos(theta_1)*sin(theta_1);
dBy1 = mag / (r1^3) * (3*cos(theta_1)^2 - 1);
dBx2 = mag / (r2^3) * 3*cos(theta_2)*sin(theta_2);
dBy2 = mag / (r2^3) * (3*cos(theta_2)^2 - 1);
% earth magnetism
Bex = Be*cos(theta_e);
Bey = Be*sin(theta_e);
% total magnetic magnitude
Bx1 = Bex + dBx1;
By1 = Bey + dBy1;
B_1 = [Bx1 ; By1];
Bx2 = Bex + dBx2;
By2 = Bey + dBy2;
B_2 = [Bx2 ; By2];
% trust index
Q = (B_1/Be)'*(B_2/Be);
V = Q - 1;
% orientation error
phi_e = atan2(Bey,Bex);
phi_1 = atan2(By1,Bx1);
phi_2 = atan2(By2,Bx2);
phi_avg = (phi_1 + phi_2) / 2;
phi_errd = (phi_avg - phi_e)*180/pi;
if phi_errd < -300
    phi_errd = phi_errd + 300;
end
%
Ans = [phi_errd V rm theta_m theta_d theta_e];
R(i,:) = Ans;
end
plot( R(:,1) , R(:,2) ,'*');