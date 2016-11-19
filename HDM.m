%%
clear all ; close all ; clc;

A = [1 2];
B = [3 4];
C = [5 6 7];
D = [8 9 10];

[X1, X2, X3, X4] = ndgrid(A, B, C, D);

c = 0;
for i = 1 : 2
    for j = 1 : 2
        for k = 1 : 3
            for z = 1 : 3
                [X1(i,j,k,z), X2(i,j,k,z), X3(i,j,k,z), X4(i,j,k,z)]
                c = c + 1
            end
        end
    end
end

%% one case
clear all ; close all ; clc;

N = 50;
N2 = 360;
Nt = N*N*N2*N2;
gap = 0.15;
sd = 0.34;
sdmin = 0.29;
sdmax = 0.39;

rmin = 0.45;
rmax = 2.5;
r = zeros(1,N);

for i = 1 : N
    r(1,i) = sqrt((rmax^2 - rmin^2) / N * (i - 1) + rmin^2);
end

r1 = r;
r2 = r;
theta1 = linspace(0, 2*pi, N2);
theta2 = linspace(0, 2*pi, N2);

[X1, X2, X3, X4] = ndgrid(r1, r2, theta1, theta2);
X = zeros(Nt,4);

c1 = 0;
for i = 1 : N
    for j = 1 : N
        for k = 1 : N2
            for z = 1 : N2
                c1 = c1 + 1;
                X(c1,:) = [X1(i,j,k,z), X2(i,j,k,z), X3(i,j,k,z), X4(i,j,k,z)];
            end
        end
    end
end

mag = 1.328*10^-7;                        % magnetic dipole
Be = 3.64*10^-5;                          % earth magnetism
theta_e = 2*pi*rand();
Bex = Be*cos(theta_e);
Bey = Be*sin(theta_e);
theta_d = 2*pi*rand();

dQ = zeros(Nt,1);
maxerr = zeros(180000,5);

c2 = 0;
c3 = 0;
for i = 1 : Nt
    if abs(X(i,3)-X(i,4)) < pi/2
        r12 = sqrt(X(i,1)^2 + X(i,2)^2 - 2*X(i,1)*X(i,2)*cos(X(i,3)-X(i,4)));
        
        if r12 >= sdmin && r12 <= sdmax
            c2 = c2 + 1;
            
            dBx1 = mag / (X(i,1)^3) * 3*cos(X(i,3))*sin(X(i,3));
            dBy1 = mag / (X(i,1)^3) * (3*cos(X(i,3))^2 - 1);
            dB1 = sqrt(dBx1^2 + dBy1^2);
            
            dBx2 = mag / (X(i,2)^3) * 3*cos(X(i,4))*sin(X(i,4));
            dBy2 = mag / (X(i,2)^3) * (3*cos(X(i,4))^2 - 1);
            dB2 = sqrt(dBx2^2 + dBy2^2);
            
            RR = [cos(-theta_d) -sin(-theta_d) ; sin(-theta_d) cos(-theta_d)];
            dBv1 = RR*[dBx1 ; dBy1];
            dBv2 = RR*[dBx2 ; dBy2];
            
            Bx1 = Bex + dBv1(1,1);
            By1 = Bey + dBv1(2,1);
            Bx2 = Bex + dBv2(1,1);
            By2 = Bey + dBv2(2,1);
            phi_e = atan2(Bey,Bex);
            phi_1 = atan2(By1,Bx1);
            phi_2 = atan2(By2,Bx2);
            phi_avg = (phi_1 + phi_2) / 2;
            phi_errd = (phi_avg - phi_e)*180/pi;
            
            dQ(i,1) = (Bex / (Be^2))*(dBv1(1,1) + dBv2(1,1)) + (Bey / (Be^2))*(dBv1(2,1) + dBv2(2,1));
            
            if dQ(i,1) <= gap && dQ(i,1) >= -gap
                if dB1 > 0.03*Be || dB2 > 0.03*Be
                    c3 = c3 + 1;
                    maxerr(c3,:) = [phi_e phi_1 phi_2 phi_avg phi_errd];
                end
            end
        end
    end
end

P1 = c3 / c2
max(maxerr(:,5))

%% many cases
clear all ; close all ; clc;

N = 50;
N2 = 360;
Nt = N*N*N2*N2;
gap = 0.15;
sd = 0.34;
sdmin = 0.29;
sdmax = 0.39;

rmin = 0.45;
rmax = 2.5;
r = zeros(1,N);

for i = 1 : N
    r(1,i) = sqrt((rmax^2 - rmin^2) / N * (i - 1) + rmin^2);
end

r1 = r;
r2 = r;
theta1 = linspace(0, 2*pi, N2);
theta2 = linspace(0, 2*pi, N2);

[X1, X2, X3, X4] = ndgrid(r1, r2, theta1, theta2);
X = zeros(Nt,4);

c1 = 0;
for i = 1 : N
    for j = 1 : N
        for k = 1 : N2
            for z = 1 : N2
                c1 = c1 + 1;
                X(c1,:) = [X1(i,j,k,z), X2(i,j,k,z), X3(i,j,k,z), X4(i,j,k,z)];
            end
        end
    end
end

mag = 1.328*10^-7;                        % magnetic dipole
Be = 3.64*10^-5;                          % earth magnetism
N3 = 3;

theta_e = linspace(0, 2*pi, N3);
theta_d = linspace(0, 2*pi, N3);

[X5, X6] = ndgrid(theta_e, theta_d);
X_ = zeros(N3*N3,2);

c2 = 0;
for i = 1 : N3
    for j = 1 : N3
        c2 = c2 + 1;
        X_(c2,:) = [X5(i,j), X6(i,j)];
    end
end

dQ = zeros(Nt,1);
P = zeros(c2,1);

c3 = 0;
c4 = 0;
for j = 1 : c2
    Bex = Be*cos(X_(j,1));
    Bey = Be*sin(X_(j,1));
    
    for i = 1 : Nt
        if abs(X(i,3)-X(i,4)) < pi/2
            r12 = sqrt(X(i,1)^2 + X(i,2)^2 - 2*X(i,1)*X(i,2)*cos(X(i,3)-X(i,4)));
            
            if r12 >= sdmin && r12 <= sdmax
                c3 = c3 + 1;
                
                dBx1 = mag / (X(i,1)^3) * 3*cos(X(i,3))*sin(X(i,3));
                dBy1 = mag / (X(i,1)^3) * (3*cos(X(i,3))^2 - 1);
                dB1 = sqrt(dBx1^2 + dBy1^2);
                
                dBx2 = mag / (X(i,2)^3) * 3*cos(X(i,4))*sin(X(i,4));
                dBy2 = mag / (X(i,2)^3) * (3*cos(X(i,4))^2 - 1);
                dB2 = sqrt(dBx2^2 + dBy2^2);
                
                RR = [cos(-X_(j,2)) -sin(-X_(j,2)) ; sin(-X_(j,2)) cos(-X_(j,2))];
                dBv1 = RR*[dBx1 ; dBy1];
                dBv2 = RR*[dBx2 ; dBy2];
                
                dQ(i,1) = (Bex / (Be^2))*(dBv1(1,1) + dBv2(1,1)) + (Bey / (Be^2))*(dBv1(2,1) + dBv2(2,1));
                
                if dQ(i,1) <= gap && dQ(i,1) >= -gap
                    if dB1 > 0.05*Be || dB2 > 0.05*Be
                        c4 = c4 + 1;
                    end
                end
            end
        end
    end
    P(j,1) = c4 / c3;
end
sum_P = 0;
for j = 1 : c2
    sum_P = sum_P + P(j,1)
end

P1 = sum_P / c2