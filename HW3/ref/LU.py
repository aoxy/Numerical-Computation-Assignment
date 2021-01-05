import numpy as np
import pandas as pd
import copy
np.random.seed(2)


def LU_decomposition(A):
    n = len(A[0])
    L = np.zeros((n, n))
    U = np.zeros((n, n))
    for i in range(n):
        L[i][i] = 1
        if i == 0:
            U[0][0] = A[0][0]
            for j in range(1, n):
                U[0][j] = A[0][j]
                L[j][0] = A[j][0] / U[0][0]
        else:
            for j in range(i, n):  # U
                temp = 0
                for k in range(0, i):
                    temp = temp + L[i][k] * U[k][j]
                U[i][j] = A[i][j] - temp
            for j in range(i + 1, n):  # L
                temp = 0
                for k in range(0, i):
                    temp = temp + L[j][k] * U[k][i]
                L[j][i] = (A[j][i] - temp) / U[i][i]
    return L, U


def LU_PDM_AUX2(A, L, U):
    n = len(A)  # A的维度
    A1 = A[1:, 1:]
    L1 = L[1:, 1:]
    U1 = U[1:, 1:]
    if n > 1:
        for i in range(1, n):
            for j in range(i, n):
                temp = A[i, j]-A[i, 0]*A[0, j]/A[0, 0]
                A[i, j] = temp
                A[j, i] = temp
        LU_PDM_AUX2(A1, L1, U1)


def LU_PDM_AUX(L, U):
    n = len(U)
    U1 = U[1:, 1:]
    if n > 1:
        for i in range(1, n):
            for j in range(i, n):
                temp = U[i, j]-U[i, 0]*U[0, j]/U[0, 0]
                U[i, j] = temp
                U[j, i] = temp
        # U[1:, 0] = np.zeros_like(U[1:, 0])
        LU_PDM_AUX(L, U1)


def LU_PDM(A):
    """
    正定矩阵(Positive Definite Matrix)LU分解算法
    """
    n = len(A)  # A的维度
    L = np.identity(n)  # 初始化L
    U = copy.deepcopy(A)  # 初始化U
    LU_PDM_AUX(L, U)
    return L, U


def dool(A):
    A = copy.deepcopy(A)
    n = A.shape[0]
    A[1:n, 0] = A[1:n, 0]/A[0, 0]
    for t in range(1, n-1):
        A[t, t:] = A[t, t:]-np.dot(A[t, :t], A[:t, t:])
        A[t+1:, t] = (A[t+1:, t]-np.dot(A[t+1:, :t], A[:t, t]))/A[t, t]
    A[n-1, n-1] = A[n-1, n-1]-np.dot(A[n-1, 0:n-1], A[0:n-1, n-1])
    L = np.tril(A, -1)+np.identity(n)
    U = np.triu(A)
    return A, L, U


def dool2(A):
    A = copy.deepcopy(A)
    n = A.shape[0]
    for i in range(1, n):
        A[i, 0] = A[i, 0]/A[0, 0]
    # for i in range(1,n-1):
    #     for j in range(j,n-1):
    #         A[i,j]=A[i,j]-
    for t in range(1, n-1):
        A[t, t:] = A[t, t:]-np.dot(A[t, :t], A[:t, t:])
        A[t+1:, t] = (A[t+1:, t]-np.dot(A[t+1:, :t], A[:t, t]))/A[t, t]
    print("dool2")
    print(A)
    A[n-1, n-1] = A[n-1, n-1]-np.dot(A[n-1, 0:n-1], A[0:n-1, n-1])
    L = np.tril(A, -1)+np.identity(n)
    U = np.triu(A)
    return A, L, U


def Doolittle(A):
    n = A.shape[0]
    L = np.zeros_like(A)
    U = np.zeros_like(A)
    for k in range(n):
        pass


if __name__ == '__main__':
    # A = np.random.randint(1, 10, size=[3, 3])  # 注意A的顺序主子式大于零
    # A = [[1, 2, 3, 4], [3, 4, 1, 2], [4, 1, 2, 3], [2, 3, 4, 1]]
    # A = [[1, -2, -2, -3], [3, -9, 0, -9], [-1, 2, 4, 7], [-3, -6, 26, 2]]
    A = [[4, -2, 4, 2], [-2, 10, -2, -7], [4, -2, 8, 4], [2, -7, 4, 7]]
    # https://blog.csdn.net/ibunny/article/details/79414232 的例子
    # A = [[1, 2, 3],
    #      [2, 3, 4],
    #      [3, 4, 5]]
    L, U = LU_decomposition(A)
    print("L矩阵是：\n", L)
    print("U矩阵是：\n", U)
    print(L.dot(U))
    B = np.array([[1, 2, 3],
                  [2, 3, 4],
                  [3, 4, 5]])
    C = np.array([[1, -2, -2, -3],
                  [3, -9, 0, -9],
                  [-1, 2, 4, 7],
                  [-3, -6, 26, 2]])
    D = np.array([[4, -2, 4, 2],
                  [-2, 10, -2, -7],
                  [4, -2, 8, 4],
                  [2, -7, 4, 7]])
    L, U = LU_PDM(D)
    print("my")
    print(U)
    print(L)
    print("Doo")
    a, b, c = dool(D)
    a2, b2, c2 = dool2(D)
    if (a == a2).all():
        print("a==a2")
    else:
        print("a=\n", a)
        print("a2=\n", a2)
    if (b == b2).all():
        print("b==b2")
    else:
        print("b=\n", b)
        print("b2=\n", b2)
    if (c == c2).all():
        print("c==c2")
    else:
        print("c=\n", c)
        print("c2=\n", c2)
    a, b, c = dool(D)
    print(b)
    print(c)

# A =
#      1    -2    -2    -3
#      3    -3     6     0
#     -1     0     2     4
#     -3     4    -2     1
# L =
#      1     0     0     0
#      3     1     0     0
#     -1     0     1     0
#     -3     4    -2     1
# U =
#      1    -2    -2    -3
#      0    -3     6     0
#      0     0     2     4
#      0     0     0     1


# [[ 1 -2 -2 -3]
#  [ 3 -3  6  0]
#  [-1  0  2  4]
#  [-3  4 -2  1]]
# [[ 1.  0.  0.  0.]
#  [ 3.  1.  0.  0.]
#  [-1.  0.  1.  0.]
#  [-3.  4. -2.  1.]]
# [[ 1 -2 -2 -3]
#  [ 0 -3  6  0]
#  [ 0  0  2  4]
#  [ 0  0  0  1]]

# def tmq(A):
#     n=len(A)
#     for i in range(n):   #计算A^(i)
#         for ix in range(i+1,n):
#             for in range(ix,n):
#         for iy=ix :n
#             A(ix)(iy)=A(ix)(iy)-A(i)(iy)*A(ix)(i)/A(i)(i)
#             A(iy)(ix)=A(ix)(iy)     %对称性
#         end
#     end
#     for k=i :n
#         u(i)(k)=A(i)(k)
#         l(k)(i)=A(i)(k)/u(i)(i)
#     #计算出 L,U。
# end

# def test_tmq():
#     print("=============tmq:=============")

# test_tmq()
