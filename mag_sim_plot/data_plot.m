%% data plot
clear all ; close all ; clc;

%%
clear all ; close all ; clc;

vx = load('vx.txt');
ut_05 = load('ut_05.txt');
st_05 = load('st_05.txt');
ut_10 = load('ut_10.txt');
st_10 = load('st_10.txt');
ut_20 = load('ut_20.txt');
st_20 = load('st_20.txt');
ut_40 = load('ut_40.txt');
st_40 = load('st_40.txt');

figure;
plot(vx,ut_05,'LineWidth',2,'Color','m');
grid on;
xlabel('Q') ; ylabel('orientation error mean (deg)');
hold on;
plot(vx,ut_10,'LineWidth',2,'Color','b');
plot(vx,ut_20,'LineWidth',2,'Color','g');
plot(vx,ut_40,'LineWidth',2,'Color','r');
legend('0.5*M_{0}','M_{0}','2*M_{0}','4*M_{0}');

figure;
plot(vx,st_05,'LineWidth',2,'Color','m');
grid on;
xlabel('Q') ; ylabel('orientation error variance (deg)');
hold on;
plot(vx,st_10,'LineWidth',2,'Color','b');
plot(vx,st_20,'LineWidth',2,'Color','g');
plot(vx,st_40,'LineWidth',2,'Color','r');
legend('0.5*M_{0}','M_{0}','2*M_{0}','4*M_{0}');
