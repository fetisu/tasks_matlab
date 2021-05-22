function Q = ElectroStaticBalls(XYZ,R,F)
XYZ = XYZ';
% Проверка неперекрытия шаров
N = length(R);
for ii = 1:(N - 1)
    for jj = (ii + 1):N
        if norm(XYZ(:,ii)-XYZ(:,jj)) <= R(ii) + R(jj)
            error('Impossible data')
        end
    end
end
% Составляем матрицу для СЛУ
M = zeros(N);
for ii = 1:N
    M(ii, ii) = 1 / R(ii);
end
for ii = 1:N
    for jj = 1:N
        if(ii ~= jj)
            M(ii, jj) = 1 / norm(XYZ(:,ii)-XYZ(:,jj));
        end
    end
end
            
Q = M \ F;
