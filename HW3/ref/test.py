import numpy as np

A = np.array([[4, -2, 4, 2],
              [-2, 10, -2, -7],
              [4, -2, 8, 4],
              [2, -7, 4, 7]])
L = np.array([[2,   0, 0,   0],
              [-1,  3,  0,    0],
              [2,  0,  2,    0],
              [1,  -2,    1,  1]])

Lt=np.transpose(L)
R=np.dot(L,Lt)
print(R)