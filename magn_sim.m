%%
clear all ; close all ; clc;

%%  plot B vector
clear all ; close all ; clc;

% given from random determine
rm = 0.5;                         % magnetic dipole center to sensor distance
theta_m = pi/3;                   % magnetic dipole center to sensor orientation
theta_d = pi/3;                   % magnetic dipole orientation
theta_e = pi/2;                   % earth magnetism orientation
% rm = 1.2171;
% theta_m = 4.5339;
% theta_d = 1.9483;
% theta_e = 3.1346;
% given from measurement
mag = 0.005;                      % magnetic dipole
Be = 0.25;                        % earth magnetism
% sensor coordinates
s1 = [0 ; 0];
s2 = [-0.11 ; -0.34];
% magnetic dipole coordinates
xm = rm*cos(theta_m);
ym = rm*sin(theta_m) - 0.14;      % from car center to magnetic dipole
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
if phi_e < 0
    phi_e = phi_e + 2*pi;
end
phi_1 = atan2(By1,Bx1);
if phi_1 < 0
    phi_1 = phi_1 + 2*pi;
end
phi_2 = atan2(By2,Bx2);
if phi_2 < 0
    phi_2 = phi_2 + 2*pi;
end
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
plot(m2(1,1),m2(2,1),'-o','MarkerSize',5,'color','r');
plot([m(1,1) m2(1,1)],[m(2,1) m2(2,1)],'LineWidth',1,'Color','r');
% two sensor to magnetic dipole line
plot([m(1,1) s1(1,1)],[m(2,1) s1(2,1)],'LineWidth',1.3,'Color','m');
plot([m(1,1) s2(1,1)],[m(2,1) s2(2,1)],'LineWidth',1.3,'Color','m');
% delta magnetic dipole vector
plot([s1(1,1) s1(1,1)+dBx1],[s1(2,1) s1(2,1)+dBy1],'LineWidth',2,'Color','k');
plot([s2(1,1) s2(1,1)+dBx2],[s2(2,1) s2(2,1)+dBy2],'LineWidth',2,'Color','k');
% earth magnetism vector
plot([s1(1,1) s1(1,1)+Bex],[s1(2,1) s1(2,1)+Bey],'-.','LineWidth',1.5,'Color','c');
plot([s2(1,1) s2(1,1)+Bex],[s2(2,1) s2(2,1)+Bey],'-.','LineWidth',1.5,'Color','c');
% total magnetic dipole vector
plot([s1(1,1) s1(1,1)+Bx1],[s1(2,1) s1(2,1)+By1],'LineWidth',1.5,'Color','g');
plot([s2(1,1) s2(1,1)+Bx2],[s2(2,1) s2(2,1)+By2],'LineWidth',1.5,'Color','g');

%%  random simulation
clear all ; close all ; clc;

% start random determine
N = 20000000;
R = zeros(N,8);

div = 100;
rV = 0.3;                    % range of V : -0.3 ~ 0.3
gap = rV/div;
n = zeros(2*div,1);
ut = zeros(2*div,1);
st = zeros(2*div,1);

for i = 1 : N
    while 1
        rxm = 3*rand();
        rym = 3*rand();
        if rxm <= 0.6 && rym <= 0.6
            ;
        else
            break;
        end
    end
    % given from random determine
    rm = sqrt( rxm^2 + rym^2 );
    theta_m = 2*pi*rand();
    theta_d = 2*pi*rand();
    theta_e = 2*pi*rand();
    % given from measurement
    mag = 0.005;                       % magnetic dipole
    Be = 0.25;                         % earth magnetism
    % sensor coordinates
    s1 = [0 ; 0];
    s2 = [-0.11 ; -0.34];
    % magnetic dipole coordinates
    xm = rm*cos(theta_m);
    ym = rm*sin(theta_m) - 0.14;       % from car center to magnetic dipole
    m = [xm ; ym];
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
    if phi_errd > 180
        phi_errd = phi_errd - 360;
    end
    if phi_errd < -180
        phi_errd = phi_errd + 360;
    end
    if phi_errd > 150
        phi_errd = phi_errd - 180;
    end
    if phi_errd < -150
        phi_errd = phi_errd + 180;
    end
    
    Ans = [phi_errd V rxm rym rm theta_m theta_d theta_e];
    R(i,:) = Ans;
    
    if V > -rV && V < rV
        k = ceil(V/gap);
        k = k + div;
        n(k,1) = n(k,1) + 1;
        ut(k,1) = ut(k,1) + phi_errd;
        st(k,1) = st(k,1) + phi_errd^2;
    end
end

for j = 1 : 2*div
    ut(j,1) = ut(j,1) / n(j,1);
    st(j,1) = st(j,1) / n(j,1);
end

plot(R(:,1),R(:,2),'*');
grid on;
xlabel('error (deg)') ; ylabel('trust index');

% figure;
% plot(R(:,3),R(:,4),'o');
% xlabel('xm (m)') ; ylabel('ym (m)');

Vx = -rV+gap : gap : rV;
figure;
plot(Vx,st);
grid on;
axis([-0.32 0.32 0 46]);
xlabel('trust index range') ; ylabel('variance');

figure;
plot(Vx,ut);
grid on;
axis([-0.32 0.32 -0.5 0.5]);
xlabel('trust index range') ; ylabel('average');

