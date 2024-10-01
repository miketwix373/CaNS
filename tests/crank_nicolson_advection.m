function u_final = crank_nicolson_advection(u0, U, dx, dt, nt)
% Solves the advection equation using the Crank-Nicolson method.
%
% Inputs:
%   u0: Initial condition as a vector.
%   U: Advection velocity.
%   dx: Spatial step size.
%   dt: Temporal step size.
%   nt: Number of time steps.
%
% Outputs:
%   u_final: Solution at the final time step.

nx = length(u0);
u = zeros(nx, nt+1);
u(:, 1) = u0;

alpha = U * dt / (4 * dx);

for j = 1:nt
    A = diag(ones(nx, 1) + alpha) - diag(alpha(2:end), 1) - diag(alpha(1:end-1), -1);
    b = zeros(nx, 1);
    b(2:end-1) = u(2:end-1, j) - alpha .* (u(3:end, j) - u(1:end-2, j));

    u(:, j+1) = A \ b;
end

u_final = u(:, end);
end