function [F,X,Y,P] = SphereDipPotential(XYZ,Q,D, R,r0,a,b,Dx,Dy,Nxy)
XYZ = XYZ';
% ВАЖНО делаем из b по оси OY
b = b - a .* dot(a,b) ./ dot(a,a);

% Матрица перехода
P = [a(1) b(1); a(2) b(2); a(3) b(3)];
dx = (Dx(2) - Dx(1)) / (Nxy(2) - 1);
dy = (Dy(2) - Dy(1)) / (Nxy(1) - 1);
N = length(R);

% Узлы для расчетов
X = [Dx(1):dx:Dx(2)];
Y = [Dy(1):dy:Dy(2)];

% Расчет потенциала в узлах
F = zeros(Nxy(1),Nxy(2));
for ii = 1:Nxy(2)
    for jj = 1:Nxy(1)
        x0 = r0(1) + a(1) * X(ii) + b(1) * Y(jj);
        y0 = r0(2) + a(2) * X(ii) + b(2) * Y(jj);
        z0 = r0(3) + a(3) * X(ii) + b(3) * Y(jj);
        r_original = [x0; y0; z0];
        for k = 1:N
            rk = r_original - XYZ(:, k);
            if(norm(rk) > R(k))
                F(jj,ii) = F(jj,ii) + Q(k) / norm(rk);
                F(jj,ii) = F(jj,ii) + dot(rk, D(k,:)') / norm(rk)^3;
            else
                F(jj,ii) = F(jj, ii) + Q(k) / R(k) + dot(rk, D(k,:)') / R(k)^3;
            end
        end
    end
end      

% Создание сетки
X = repmat(X, Nxy(1), 1);
Y = Y';
Y = repmat(Y, 1, Nxy(2));
end


