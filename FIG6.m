clc
clear
t=0:0.1:20000; 

r = 4.81; d = 0.0289; m = 2; a1 = 0.0395; b1 = 0.6; bc = 0.0451;
ka = 1.316; kb = 0.045; pr = 1.77; pc = 1.23; n1 = 2; n2 = 1; Bn = 1.5; Ic = 1;
y0 = [4.3; 0.1]; 


figure;

% tau = 50
subplot(2, 1, 1)
tau = 50; 
f = @(t, y, Z) [r*y(2)-d*y(1);
    ((a1 + b1) * ((ka * Z(1) + 1)^n1 + (pr / d)^n1) / ((ka * Z(1) + 1)^n1 * (kb * Bn + 1) + (pr / d)^n1) - b1 - bc * (Ic * (pc / d)^n2) / (Ic * (pc / d)^n2 + (ka * y(1) + 1)^n2)) * y(2) * (1 - y(2)^m)
];
sol = dde23(f, tau, y0, t);

% Plotting a 3D phase diagram with the z-axis at time
plot3(sol.x, sol.y(1, :), sol.y(2, :), 'LineWidth', 1.5);
grid on;
box on;
set(gca, 'LineWidth', 2, 'FontSize', 13);
set(gca, 'GridLineStyle', ':', 'LineWidth', 1);
xlabel('$t$', 'Interpreter', 'LaTex', 'FontSize', 13);
ylabel('$A$', 'Interpreter', 'LaTex', 'FontSize', 13);
zlabel('$N$', 'Interpreter', 'LaTex', 'FontSize', 13);
title('3D phase plot with $\tau = 50$', 'Interpreter', 'LaTex', 'FontSize', 13);
text(0.98, 0.7, '(a)', 'Units', 'normalized', 'FontSize', 13, 'HorizontalAlignment', 'right');
hold on;

% tau = 60
subplot(2, 1, 2)
tau = 60; 
f = @(t, y, Z) [r*y(2)-d*y(1);
    ((a1 + b1) * ((ka * Z(1) + 1)^n1 + (pr / d)^n1) / ((ka * Z(1) + 1)^n1 * (kb * Bn + 1) + (pr / d)^n1) - b1 - bc * (Ic * (pc / d)^n2) / (Ic * (pc / d)^n2 + (ka * y(1) + 1)^n2)) * y(2) * (1 - y(2)^m)
];
sol = dde23(f, tau, y0, t);


plot3(sol.x, sol.y(1, :), sol.y(2, :), 'LineWidth', 1.5);
grid on;
box on;
set(gca, 'LineWidth', 2, 'FontSize', 13);
set(gca, 'GridLineStyle', ':', 'LineWidth', 1);
xlabel('$t$', 'Interpreter', 'LaTex', 'FontSize', 13);
ylabel('$A$', 'Interpreter', 'LaTex', 'FontSize', 13);
zlabel('$N$', 'Interpreter', 'LaTex', 'FontSize', 13);
title('3D phase plot with $\tau = 60$', 'Interpreter', 'LaTex', 'FontSize', 13);
text(0.98, 0.7, '(b)', 'Units', 'normalized', 'FontSize', 13, 'HorizontalAlignment', 'right');
hold on;
