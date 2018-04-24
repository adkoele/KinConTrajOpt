function [f,dfdx,dfdu,dfdl] = dynamics(x,u,l,params)
%this is also a problem specific function ??

L = params.l;
m = params.m;
g = params.g;

f(1) = -1/m*(l*x(1)/L+x(2)/L^2*u);
f(2) = -1/m*(l*x(2)/L-x(1)/L^2*u)-g;
f = f';

dfdx = [-1/m*l/L -1/m*u/L^2;
    1/m*u/L^2 -1/m*l/L];
dfdu = [-1/m*x(2)/L^2;1/m*x(1)/L^2];
dfdl = [-1/m*x(1)/L;-1/m*x(2)/L];