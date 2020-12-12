clear, clc, clf
syms y(x)
F=dsolve(diff(y)==x*exp(-5*x)-5*y,y(0)==0,x)
fun = matlabFunction(F)