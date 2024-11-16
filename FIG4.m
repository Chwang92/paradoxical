clear;
clc;
tic;  

% Parameter Definition
r = 4.81; d = 0.0289; m = 2; a1 = 0.0395; b1 = 0.6; bc = 0.0451;
ka = 1.316; kb = 0.045; pr = 1.77; pc = 1.23; n1 = 2; n2 = 1; Bn = 1.5; Ic = 1;

% Initializing the grid
scale = 100;
[x, y] = meshgrid(0:0.5:200, 0:0.01:1); % Increase grid spacing to reduce flow density

u = zeros(size(x));
v = zeros(size(y));

% Calculate the direction of each point on the grid
for k = 1:numel(x)
    [F, ~] = Fdydx(0, [x(k); y(k)], r, d, m, a1, b1, bc, ka, kb, pr, pc, n1, n2, Bn, Ic);
    u(k) = F(1);
    v(k) = F(2);
end

% Mapping of flow lines
figure();
streamslice(x, y*scale, u, v*scale, 2);
set(gca,'LineWidth',2,'FontSize',13);
set(gca, 'GridLineStyle' ,':','linewidth',1');
xlabel('$A$','Interpreter','LaTex','FontSize',13);
ylabel('$N$','Interpreter','LaTex','FontSize',13);
box on;

% Marking the equilibrium of the system
syms xs ys
dxdt_eq = r*ys - d*xs;
dydt_eq = ((a1 + b1)*((ka*xs + 1)^n1 + (pr/d)^n1)/((ka*xs + 1)^n1*(kb*Bn + 1) + (pr/d)^n1) - b1 - bc*(Ic*(pc/d)^n2)/(Ic*(pc/d)^n2 + (ka*xs + 1)^n2))*ys*(1 - ys^m);

% Solve the equation to get the equilibrium point
eqns = [dxdt_eq == 0, dydt_eq == 0];
S = solve(eqns, [xs, ys], 'Real', true);

% Define different colors for equilibrium
colors = ['r', 'g', 'g', 'r','y', 'c']; 
color_index = 1;

% Mark the equilibrium point 
hold on;
for i = 1:length(S.xs)
    x_val = double(S.xs(i));
    y_val = double(S.ys(i)) * scale;  % Apply the same y-scaling
    if x_val >= 0 && x_val <= 200 && y_val >= 0 && y_val <= 1*scale
        % Select color cycling
        plot_color = colors(mod(color_index + 1, length(colors)) - 1);
        plot(x_val, y_val, 'o', 'MarkerSize', 10, 'MarkerFaceColor', plot_color);
        
        % Add Tags
        text(x_val, y_val, sprintf('Point % d', i), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', plot_color);
        
        % Update color index
        color_index = color_index +  1;
    end
end
hold off;

 yticks([0*scale 0.25*scale 0.5*scale 0.75*scale 1*scale])
 yticklabels(num2str(yticks'/scale))
xlim([0 200])
ylim([0*scale 1*scale])
fprintf('Total computation time: % .2f seconds.\n', toc);

% Define the equations of the system
function [F, Output] = Fdydx(~, XY, r, d, m, a1, b1, bc, ka, kb, pr, pc, n1, n2, Bn, Ic)
    % Extracting Variables
   x = XY(1);
   y = XY(2);

    % Calculating dx/dt and dy/dt
    dxdt = r*y - d*x;
    dydt = ((a1 + b1)*((ka*x + 1)^n1 + (pr/d)^n1)/((ka*x + 1)^n1*(kb*Bn + 1) + (pr/d)^n1) - b1 - bc*(Ic*(pc/d)^n2)/(Ic*(pc/d)^n2 + (ka*x + 1)^n2))*y*(1 - y^m);

    % Return results
    F = [dxdt; dydt];
    Output = [];
end