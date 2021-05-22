function [F,X,Y,P] = SpherePotential(XYZ,Q,R,r0,a,b,Dx,Dy,Nxy)
XYZ = XYZ';
% ВАЖНО делаем из b по оси OY
b = b - a .* dot(a,b) ./ dot(a,a);

% Матрица перехода
P = zeros(2,3);
P(1,:) = a;
P(2,:) = b;
dx = (Dx(2) - Dx(1)) / (Nxy(2) - 1);
dy = (Dy(2) - Dy(1)) / (Nxy(1) - 1);


% Узлы для расчетов
X = Dx(1):dx:Dx(2);
Y = Dy(1):dy:Dy(2);

% Расчет потенциала в узлах
F = zeros(Nxy(1),Nxy(2));
for i = 1:Nxy(2)
    for j = 1:Nxy(1)
        x0 = r0(1) + a(1) * X(i) + b(1) * Y(j);
        y0 = r0(2) + a(2) * X(i) + b(2) * Y(j);
        z0 = r0(3) + a(3) * X(i) + b(3) * Y(j);
        r_original = [x0; y0; z0];
        for k = 1:N
            if(norm(XYZ(:,k)-r_original) > R(k))
                F(j,i) = F(j, i) + Q(k) / norm(XYZ(:,k)-r_original);
            else
                F(j,i) = F(j, i) + Q(k) / R(k);
            end
        end
    end
end      

% Создание сетки
X = repmat(X, Nxy(1), 1);
Y = Y';
Y = repmat(Y, 1, Nxy(2));
end

