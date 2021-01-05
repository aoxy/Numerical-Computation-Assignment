%% 正定矩阵(Positive Definite Matrix)LU分解算法
function [L, U] = LU_PDM(A)
    [~, n] = size(A);
    L = zeros(n, n);
    U = zeros(n, n);
    for i = 1:n
        for j = i + 1:n
            for k = j:n
                A(j, k) = A(j,k)-A(i,k)*A(j,i)/A(i,i);
                A(k, j) = A(j,k);
            end
        end
        for k = i:n
            U(i, k) = A(i, k);
            L(k, i) = A(i, k)/U(i, i);
        end
    end
end