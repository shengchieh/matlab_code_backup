clear all;
clc;

N = 50;
x = linspace(5,-5,N);
y = x;
[X,Y] = meshgrid(x);
F1 = 0*X;
F2 = 0*X;

Be = 1;
THEe = pi/2;
Bex = Be*sin(THEe);
Bey = Be*cos(THEe);
A = 1;

for i = 1:N
    for j = 1:N
        R = X(i,j)*X(i,j) +  Y(i,j)*Y(i,j);
        THE = atan2(X(i,j),Y(i,j));
        
        Bx = Bex + A*3*sin(THE)*cos(THE)/(R^3);
        By = Bey + A*(3*cos(THE)*cos(THE) - 1)/(R^3);
        
        F1(i,j) = atan2(Bx,By);
        F2(i,j) = sqrt(Bx^2+ By^2);
        if F2(i,j) > 2
            F2(i,j) = 2;
        end
    end
end

Plate = 0*F1;

figure(1)
surf(X,Y,F1);
hold on
surf(X,Y,Plate);

figure(2)
surf(X,Y,F2);
hold on
surf(X,Y,Plate);

hold off