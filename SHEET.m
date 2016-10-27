clear all; close all ; clc;

N = 100;
x = linspace(5,-5,N);
y = x;
[X,Y] = meshgrid(x,y);
F1 = 0*X;
F2 = 0*X;

Be = 1;
THEe = 0;
Bex = Be*cos(THEe);
Bey = Be*sin(THEe);
A = 1;

for i = 1:N
    for j = 1:N
        R = sqrt(X(i,j)*X(i,j) +  Y(i,j)*Y(i,j));
        THE = atan2(X(i,j),Y(i,j));
        
        Bx = Bex + A*3*cos(THE)*sin(THE) / (R^3);
        By = Bey + A*(3*cos(THE)*cos(THE) - 1) / (R^3);
        
        F1(i,j) = atan2(By,Bx);
        F2(i,j) = sqrt(Bx^2 + By^2);
        if F2(i,j) > 2
            F2(i,j) = 2;
        end
    end
end

Plate = 0*F1;
Plate = Plate + THEe;

figure(1);
surfc(X,Y,F1);
hold on;
surfc(X,Y,Plate);

figure(2);
surfc(X,Y,F2);
