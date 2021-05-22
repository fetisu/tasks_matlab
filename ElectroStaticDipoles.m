function [Q,D] = ElectroStaticDipoles(XYZ,R,F)
XYZ = XYZ';
N = length(R);
% F = F1 ... Fn, E1x ... Enx, E1y ... Eny, E1z ... Enz
F = [F; zeros(3 * N, 1)];
% ans = Q1 ... Qn, p1x ... pnx, p1y ... pny, p1z ... pnz
% Проверка неперекрытия шаров
for ii = 1:(N - 1)
    for jj = (ii + 1):N
        if norm(XYZ(:,ii)-XYZ(:,jj)) <= R(ii) + R(jj)
            error('Impossible data')
        end
    end
end

% Составляем матрицу для СЛУ
M = zeros(4 * N);

% Условия на потенциал
for ii = 1 : N
    for jj = 1 : 4 * N
        j = mod(jj, N);
        if(j == 0)
            j = N;
        end
        coord = (jj - j) / N;
        r = XYZ(1:3,ii)-XYZ(1:3,j);
        if(ii == jj)
            M(ii, jj) = 1 / R(ii);
        else
            if(coord == 0)
                M(ii, jj) = 1 / norm(r);
            else
                if(j == ii)
                    M(ii, jj) = 0;
                else
                    M(ii, jj) = r(coord) / norm(r) ^ 3;
                end
            end
        end   
    end
end

%Условия на поле в центрах

for ii = (N + 1) : (4 * N)
    for jj = 1 : (4 * N)
        j = mod(jj, N);
        i = mod(ii, N);
        if(j == 0)
            j = N;
        end
        if(i == 0)
            i = N;
        end
        coordp = (jj - j) / N;
        coordE = (ii - i) / N;
        r = XYZ(1:3,i)-XYZ(1:3,j);
        if(i == jj)
            M(ii, jj) = 0;
        else
            if(ii == jj)
                M(ii, jj) = 1 / R(i) ^ 3;
            else
                if(coordp == 0)
                    M(ii, jj) = r(coordE) / norm(r) ^ 3;
                else
                    if(j == i)
                        M(ii, jj) = 0;
                    else
                        if(coordE == coordp)
                            M(ii, jj) = (3 * r(coordE) ^ 2 - norm(r)^2 ) / norm(r) ^ 5;
                        else
                            M(ii, jj) = 3 * r(coordE) * r(coordp)/ norm(r) ^ 5;
                        end
                    end
                end
            end
        end 
    end
end
Help = M \ F;
Q = Help(1:N, 1);
Dx = Help(N+1:2*N,1)';
Dy = Help(2*N+1:3*N,1)';
Dz = Help(3*N+1:4*N,1)';
D = [Dx; Dy; Dz]';
D = -D;
end



