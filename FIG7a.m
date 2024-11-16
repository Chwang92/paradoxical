
clear all
clc

%%% Parameters
r = 4.81;
d = 0.0289;
m = 2;
a1 = 0.0395;
b1 = 0.6;
bc = 0.0451;
ka = 1.316;
kb = 0.045;
pr = 1.77;
pc = 1.23;
n1 = 2;
n2 = 1;
Bn = 1.5;
Ic = 1;

%%% Define DDE function
f = @(t, y, Z) [
    r * y(2) - d * y(1);
    ((a1 + b1) * ((ka * Z(1) + 1)^n1 + (pr / d)^n1) / ((ka * Z(1) + 1)^n1 * (kb * Bn + 1) + (pr / d)^n1) - b1 - bc * (Ic * (pc / d)^n2) / (Ic * (pc / d)^n2 + (ka * y(1) + 1)^n2)) * y(2) * (1 - y(2)^m)
];

%%% Initial conditions
y0 = [4.3; 0.1];

%%% History function (constant initial condition)
hist = @(t) y0;

%%% Time delay and time span
tau_range = 50:0.1:110;  % Delay values ranging from 50 to 110 with 0.1 increments
Tstart = 0;               % Start time of solution
Tend = 50000;              % End time of solution

%%% Preallocate max_points and min_points
max_points = nan(1, length(tau_range));
min_points = nan(1, length(tau_range));


%%% Solve DDE for different values of tau and capture dynamics
for i = 1:length(tau_range)
    
    % Set delay
    tau = tau_range(i);
    
    % Solve DDE using dde23
    opts = ddeset('RelTol', 1e-6, 'AbsTol', 1e-6);
    sol = dde23(f, tau, hist, [Tstart Tend], opts);
    
    % Capture max and min betwen t=48000 and Tend
    time_indices = find(sol.x >= 48000 & sol.x <= 50000);
        max_points(i) = max(sol.y(2, time_indices));
        min_points(i) = min(sol.y(2, time_indices));
end

%%% Plot bifurcation diagram
figure(1),
plot(tau_range, max_points, 'b', 'LineWidth', 2), hold on
plot(tau_range, min_points, 'b', 'LineWidth', 2), hold on
xlabel('$\tau$','Interpreter','LaTex')
ylabel('$N$','Interpreter','LaTex')
plot(58.95, 0.28, 's', 'MarkerSize', 12, 'MarkerFaceColor', 'k'), hold on
plot(99, 0, 'd', 'MarkerSize', 12, 'MarkerFaceColor', 'k'), hold on
plot(tau_range(tau_range >= 58.95 & tau_range < 99), 0.28 * ones(1, length(tau_range(tau_range >= 58.95 & tau_range < 99))), 'k--', 'LineWidth', 2), hold on
yticks([0 0.2 0.4 0.6 0.8 1])
xticks([50 60 70 80 90 100 110])
set(gca, 'FontSize', 20);

%------------------------------------------------------------------------

