%%
clear all ; close all ; clc;

%%
clear all ; close all ; clc;

r = 0.078;        % ½ü¥b®|
w = 0.393;        % ¨â½ü¶¡¶Z
n = 1834;         % encoder count per revolution = 28 , gear ratio = 65.5 , 28*65.5=1834

%% 
clear all ; close all ; clc;

A = [0 1 ; 0 0];
B = [1 ; 0];
C = [1 0];
R = 5;
theta = 3*pi/4;
p = [R*cos(theta)+i*sin(theta) R*cos(theta)-i*sin(theta)];
L = place(A', C', p)';
check = eig(A-L*C)
ALC = A-L*C;

%%
figure(1);
plot(tout,o_t,'LineWidth',1.5,'Color','b');
grid on;
xlabel('t (sec)') ; ylabel('true omega (rad/sec)');

figure(2);
plot(tout,t_t,'LineWidth',1.5,'Color','b');
grid on;
xlabel('t (sec)') ; ylabel('true theta (rad)');

figure(3);
plot(tout,y,'LineWidth',0.1,'Color','b');
grid on;
xlabel('t (sec)') ; ylabel('y (rad)');

figure(4);
plot(tout,o_b,'LineWidth',0.1,'Color','b');
grid on;
xlabel('t (sec)') ; ylabel('omega drift (rad/sec)');

figure(5);
plot(tout,o_g,'LineWidth',0.1,'Color','b');
grid on;
xlabel('t (sec)') ; ylabel('omega gyro (rad/sec)');

figure(6);
plot(tout,o_b,'LineWidth',0.1,'Color','b');
grid on ; hold on;
plot(tout,b_h,'LineWidth',0.1,'Color','r');
xlabel('t (sec)') ; ylabel('estimated drift (rad/sec)');
legend('omega drift','estimated drift');

figure(7);
plot(tout,t_t,'LineWidth',3.5,'Color','b');
grid on ; hold on;
plot(tout,t_h,'LineWidth',1.5,'Color','r');
xlabel('t (sec)') ; ylabel('estimated theta (rad)');
legend('true theta','estimated theta');
