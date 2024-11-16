clc;
clear;
t = 0:3000;
% Setting parameters
r = 4.81; d = 0.0289; m = 2; a1 = 0.0395; b1 = 0.6; bc = 0.0451;
ka = 1.316; kb = 0.045; pr = 1.77; pc = 1.23; n1 = 2; n2 = 1; Bn = 1.5; Ic = 1; 

% Define two sets of initial conditions
y0s = {[4.3; 0.01], [4.3; 0.1]};

% Creating a Graphics Window
figure;

% Create an empty array to store tau values for color bar
tau_values = [];

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
        
        % Store tau for the color bar
        tau_values = [tau_values, tau];
        
        % Setting the color based on the value of tau
        R = sin(pi * tau / 320);    
        G = sin(pi * tau / 320 + pi/2); 
        B = cos(pi * tau / 320);   

        R = max(0, min(1, R));
        G = max(0, min(1, G));
        B = max(0, min(1, B));

        % plot the solution with the color based on tau
        plot(sol.x, sol.y(2, :), 'LineWidth', 2, 'Color', [R, G, B]);
        hold on;
    end

    % Add a color bar to represent tau values
    colormap jet; % Choose a colormap (jet gives a smooth transition from blue to red)
    colorbar;  % Add colorbar to show the color scale
    caxis([1 160]);  % Set the limits of the color scale (based on the range of tau)

    % Add labels and other formatting
    grid on;
    box on;
    set(gca, 'LineWidth', 2, 'FontSize', 13);
    set(gca, 'GridLineStyle', ':', 'LineWidth', 1);
    xlabel('$t$', 'Interpreter', 'LaTex', 'FontSize', 13);
    ylabel('$N$', 'Interpreter', 'LaTex', 'FontSize', 13);
    
    % Add (a) or (b) at the top-right corner of the subplot
    text(0.95, 0.9, ['(', char(96+i), ')'], 'Units', 'normalized', 'FontSize', 15, 'HorizontalAlignment', 'center');
    
    hold off;
end
