clc; clear;
root = pwd;
path = root + "/scalability/";

filename = 'sscaling_wtime.out';

data = readtable(path + filename,'FileType','text');

figure(1)
semilogx(data.Var1,data.Var2,'s','MarkerSize',10,'MarkerEdgeColor','black')
fontsize(gca, 13,'points')
title('Runtime Performance','Interpreter','latex','FontSize',25)
xlabel('$N$','Interpreter','latex','FontSize',25)
ylabel('$t_w$','Interpreter','latex','FontSize',25)
grid on
ylim([0.1,40])
xlim([100,20000])
