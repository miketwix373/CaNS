clc; clear;
root = pwd;
path = root + "/run/data/";


data = readtable(path+'prueba.csv');
x = reshape(data.Points_0,[256,128])';
z = reshape(data.Points_2,[256,128])';
v = sqrt(reshape(data.Velocity_0,[256,128])'.^2+reshape(data.Velocity_1,[256,128])'.^2+reshape(data.Velocity_2,[256,128])'.^2);

figure(1)
contourf(x,z,v)
colorbar

figure(2)
hold on
posX = floor(linspace(1,256,10));
for i=1:5
  plot(v(:,posX(i)),z(:,posX(i)))
  
end
hold off
legend(string(x(1,posX)))
fontsize(gca, 13,'points')
title('Bl profile','Interpreter','latex','FontSize',25)
xlabel('$V [m/s]$','Interpreter','latex','FontSize',25)
ylabel('$Z [m]$','Interpreter','latex','FontSize',25)
grid on

