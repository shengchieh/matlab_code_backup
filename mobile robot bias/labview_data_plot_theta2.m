%% LabVIEW data plot
clear all ; close all ; clc;

%%
clear all ; close all ; clc;

data = load('exp.txt');
t = 0.0046 : 0.0046 : length(data(:,1))*0.0046;          % send data period

figure;
plot(t,data(:,1),'r',t,data(:,2),'k',t,data(:,3),'b',t,data(:,4),'g',t,data(:,4)*0-360,'m--','linewidth',1.5);
grid on;
legend('theta gyro','theta encoder','theta magn & enc','theta fusion','-360 deg');
xlabel('sec') ; ylabel('deg');

figure;
plot(t,data(:,5)*57.3,'m',t,data(:,6)*57.3,'b','linewidth',1.5);
grid on;
legend('estimated drift','add drift');
xlabel('sec') ; ylabel('deg/sec');
hold on;
% plot(t,-data(:,6),'g--','linewidth',1);
