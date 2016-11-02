%%
clear all ; close all ; clc;

A = [1 2];
B = [3 4];
C = [5 6];
D = [7 8];

[X1, X2, X3, X4] = ndgrid(A, B, C, D);

c = 0;
for i = 1 : 2
    for j = 1 : 2
        for k = 1 : 2
            for z = 1 : 2
                [X1(i,j,k,z), X2(i,j,k,z), X3(i,j,k,z), X4(i,j,k,z)]
                c = c + 1
            end
        end
    end
end

%%
clear all ; close all ; clc;

N = 50;
Nt = N^4;
gap = 0.15;
sd = 0.33;
sdmin = 0.23;
sdmax = 0.43;

rmin = 0.45;
rmax = 2.5;
r = zeros(1,N);

for i = 1 : N
    r(1,i) = sqrt((rmax^2 - rmin^2) / N * (i - 1) + rmin^2);
end

r1 = r;
r2 = r;
theta1 = linspace(0, 2*pi, N);
theta2 = linspace(0, 2*pi, N);

[X1, X2, X3, X4] = ndgrid(r1, r2, theta1, theta2);
X = zeros(Nt,4);

c1 = 0;
for i = 1 : N
    for j = 1 : N
        for k = 1 : N
            for z = 1 : N
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

dQ = zeros(Nt,1);

c2 = 0;
c3 = 0;
for i = 1 : Nt
    if abs(X(i,3)-X(i,4)) < pi
        r12 = sqrt(X(i,1)^2 + X(i,2)^2 - 2*X(i,1)*X(i,2)*cos(X(i,3)-X(i,4)));
        
        if r12 > sdmin && r12 < sdmax
            c2 = c2 + 1;
            
            dBx1 = mag / (X(i,1)^3) * 3*cos(X(i,3))*sin(X(i,3));
            dBy1 = mag / (X(i,1)^3) * (3*cos(X(i,3))^2 - 1);
            dB1 = sqrt(dBx1^2 + dBy1^2);
            
            dBx2 = mag / (X(i,2)^3) * 3*cos(X(i,4))*sin(X(i,4));
            dBy2 = mag / (X(i,2)^3) * (3*cos(X(i,4))^2 - 1);
            dB2 = sqrt(dBx2^2 + dBy2^2);
            
            dQ(i,1) = (Bex / (Be^2))*(dBx1 + dBx2) + (Bey / (Be^2))*(dBy1 + dBy2);
            
            if dQ(i,1) <= gap && dQ(i,1) >= -gap
                if dB1 > 0.01*Be || dB2 > 0.01*Be
                    c3 = c3 + 1;
                end
            end
        end
    end
end

P1 = c3 / c2
