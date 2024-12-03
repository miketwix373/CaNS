%% Test for perturbation

clear;
clc;

x = linspace(0,2,1000);
z = linspace(0,0.5,200);
x0 = 0.1;
lx = 0.1;
ly = 0.1;

[XX,ZZ] = meshgrid(x,z);
pt = exp(-((XX-x0)./lx).^2+((ZZ)./ly).^2);

contourf(XX,flip(ZZ,1),pt)

 