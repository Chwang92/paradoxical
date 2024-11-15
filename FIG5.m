clc;
clear;
t = 0:3000;
% Setting parameters
r = 4.81; d = 0.0289; m = 2; a1 = 0.0395; b1 = 0.6; bc = 0.0451;
ka = 1.316; kb = 0.045; pr = 1.77; pc = 1.23; n1 = 2; n2 = 1; Bn = 1.5; Ic = 1; 

% Define two sets of initial conditions
y0s = {[4.3; 0.01], [4.3;0.1]};

% Creating a Graphics Window
figure;

for i = 1:length(y0s)
    y0 = y0s{i};
    subplot(2, 1, i); 
    % Loop traversal time lag value
for tau = 1:5:160
    % Defining the equation
    f = @(t, y, Z) [
        r*y(2) - d*y(1);
        ((a1 + b1)*((ka*Z(1) + 1)^n1 + (pr/d)^n1)/((ka*Z(1) + 1)^n1*(kb*Bn + 1) + (pr/d)^n1) - b1 - bc*(Ic*(pc/d)^n2)/(Ic*(pc/d)^n2 + (ka*y(1) + 1)^n2))*y(2)*(1 - y(2)^m)
       ];
    
    % solving the equation
    sol = dde23(f, tau, y0, t);
    
 % Setting the color
    R = sin(pi * tau / 320);    
    G = sin(pi * tau / 320 + pi/2); 
    B = cos(pi * tau / 320);   

    
    R = max(0, min(1, R));
    G = max(0, min(1, G));
    B = max(0, min(1, B));

    % plot
    plot(sol.x, sol.y(2, :), 'LineWidth', 2, 'Color', [R, G, B]);
    hold on;
end

grid on;
box on;
set(gca,'LineWidth',2,'FontSize',13);
set(gca, 'GridLineStyle' ,':','linewidth',1');
xlabel('Time (hours)','FontSize',13);
ylabel('$N$','Interpreter','LaTex','FontSize',13);
hold off;
end