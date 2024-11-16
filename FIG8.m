clc
clear
t = 0:0.1:5000; 

r = 4.81; d = 0.0289; m = 2; a1 = 0.0395; b1 = 0.6; bc = 0.0451;
ka = 1.316; kb = 0.045; pr = 1.77; pc = 1.23; n1 = 2; n2 = 1; Ic = 1;
y1 = [4.3; 0.03]; 
y2 = [4.3; 0.02];
tau = 72; 

figure;

% First subplot
subplot(3, 1, 1)
Bn = 1.3;
f = @(t, y, Z) [r * y(2) - d * y(1);
    ((a1 + b1) * ((ka * Z(1) + 1)^n1 + (pr / d)^n1) / ((ka * Z(1) + 1)^n1 * (kb * Bn + 1) + (pr / d)^n1) - b1 - bc * (Ic * (pc / d)^n2) / (Ic * (pc / d)^n2 + (ka * y(1) + 1)^n2)) * y(2) * (1 - y(2)^m)
];
sol = dde23(f, tau, y1, t);

plot(sol.x, sol.y(2, :), 'r','LineWidth', 1.5);
grid on;
box on;
set(gca, 'LineWidth', 2, 'FontSize', 13);
set(gca, 'GridLineStyle', ':', 'LineWidth', 1);
set(gca, 'XTick', []);
set(gca, 'XTickLabel', []);
ylabel('$N$', 'Interpreter', 'LaTex', 'FontSize', 13);
legend('$(B_n, N_0)=(1.3,0.03)$', 'Interpreter', 'LaTeX', 'Location', 'best');
hold on;
text(0.95, 0.05, '(a)', 'Units', 'normalized', 'FontSize', 14);
% Second subplot
subplot(3, 1, 2)
Bn = 1.5;
f = @(t, y, Z) [r * y(2) - d * y(1);
    ((a1 + b1) * ((ka * Z(1) + 1)^n1 + (pr / d)^n1) / ((ka * Z(1) + 1)^n1 * (kb * Bn + 1) + (pr / d)^n1) - b1 - bc * (Ic * (pc / d)^n2) / (Ic * (pc / d)^n2 + (ka * y(1) + 1)^n2)) * y(2) * (1 - y(2)^m)
];
sol = dde23(f, tau, y1, t);

plot(sol.x, sol.y(2, :), 'k','LineWidth', 1.5);
grid on;
box on;
set(gca, 'LineWidth', 2, 'FontSize', 13);
set(gca, 'GridLineStyle', ':', 'LineWidth', 1);
% Hide x-axis ticks and labels
set(gca, 'XTick', []);
set(gca, 'XTickLabel', []);
ylabel('$N$', 'Interpreter', 'LaTex', 'FontSize', 13);
legend('$(B_n,N_0)=(1.5,0.03)$', 'Interpreter', 'LaTeX', 'Location', 'best');
hold on;
text(0.95, 0.05, '(b)', 'Units', 'normalized', 'FontSize', 14);
% Third subplot
subplot(3, 1, 3)
Bn = 1.3;
f = @(t, y, Z) [r * y(2) - d * y(1);
    ((a1 + b1) * ((ka * Z(1) + 1)^n1 + (pr / d)^n1) / ((ka * Z(1) + 1)^n1 * (kb * Bn + 1) + (pr / d)^n1) - b1 - bc * (Ic * (pc / d)^n2) / (Ic * (pc / d)^n2 + (ka * y(1) + 1)^n2)) * y(2) * (1 - y(2)^m)
];
sol = dde23(f, tau, y2, t);

plot(sol.x, sol.y(2, :), 'LineWidth', 1.5);
grid on;
box on;
set(gca, 'LineWidth', 2, 'FontSize', 13);
set(gca, 'GridLineStyle', ':', 'LineWidth', 1);
xlabel('$t$', 'Interpreter', 'LaTex', 'FontSize', 13);
ylabel('$N$', 'Interpreter', 'LaTex', 'FontSize', 13);
legend('$(B_n, N_0)=(1.3,0.02)$', 'Interpreter', 'LaTeX', 'Location', 'best');
hold on;
text(0.95, 0.05, '(c)', 'Units', 'normalized', 'FontSize', 14);