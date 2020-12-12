clear, clc, clf

syms x;
syms a;
syms h;
f = (x-a-2*h)*(x-a-h)*(x-a)*(x-a+h);
int(f,x,a,a+2*h)