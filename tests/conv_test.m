%% Convective bounbdary conditions test
clc
clear
umat = [
    1, 2, 3
    2, 4, 6
    3, 6, 9
    4, 8, 12
    5, 10, 15
    6, 12, 18
    7, 14, 21
    8, 16, 24
];

for j=1:length(umat(:,1))
    U(j) = mean(umat(j,:));
end

dx = 0.1;
dt = 0.01;
u_final = advection(U,[umat(:,2),umat(:,1)], dt, dx);

