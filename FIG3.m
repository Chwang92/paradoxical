% Parameter values
gamma = 4.81;
delta = 0.0289;
m = 2;
alpha = 0.0395;
beta = 0.6;
beta_c = 0.0451;
k_a = 1.316;
k_b = 0.045;
p_r = 1.77;
p_c = 1.23;
n1 = 2;
n2 = 1;

% Defining the G(A)
G = @(A, B_n, I_c) ((alpha + beta) * ((k_a * A + 1).^n1 + (p_r / delta)^n1) ./ ...
    ((k_a * A + 1).^n1 .* (k_b * B_n + 1) + (p_r / delta)^n1) - beta - ...
    (beta_c * I_c * (p_c / delta)^n2) ./ ...
    (I_c * (p_c / delta)^n2 + (k_a * A + 1).^n2));

% Defining the derivative G'(A)
G_prime = @(A, B_n, I_c) beta_c * n2 * I_c * k_a * ((k_a * A + 1).^(n2 - 1)) .* (p_c / delta)^n2 ./ ...
    ((I_c * (p_c / delta)^n2 + (k_a * A + 1).^n2)).^2 - ...
    (alpha + beta) * n1 * k_a * k_b * B_n * ((k_a * A + 1).^(n1 - 1)) .* ((p_r / delta)^n1) ./ ...
    (((k_a * A + 1).^n1 .* (k_b * B_n + 1) + (p_r / delta)^n1)).^2;

% Define the grid range [0, 10] in steps of 0.1.
[B_n, I_c] = meshgrid(0:0.005:10, 0:0.005:10);

% Calculating G(gamma/delta) and G(0)
A1 = gamma / delta;
A2 = 0;
G1 = G(A1, B_n, I_c);  % G(gamma/delta)
G2 = G(A2, B_n, I_c);  % G(0)

% Find the region satisfying G(gamma/delta) > 0 and G(0) < 0
region1 = (G1 > 0) & (G2 < 0);

% Initialize region2 as false
region2 = false(size(B_n));

% Iteratively solve for P and check conditions
for i = 1:size(B_n, 1)
    for j = 1:size(B_n, 2)
        b_n = B_n(i, j);
        i_c = I_c(i, j);
        
        % Define a function handle for G with current B_n and I_c
        G_current = @(A) G(A, b_n, i_c);
        
        % Use an initial guess for A (e.g., gamma/delta) and search for a root of G(P) = 0
        try
            P_candidate = fzero(G_current, 0.5*gamma / delta);  % You can adjust the initial guess as needed
            
            % Filter out P that satisfies (delta * P) / gamma between 0 and 1
            if (delta * P_candidate) / gamma > 0 && (delta * P_candidate) / gamma < 1
                % Calculate G'(P)
                Gp = G_prime(P_candidate, b_n, i_c);
                
                % Check if G(0) < 0 and G'(P) < 0
                if G2(i, j) < 0 && Gp < 0
                    region2(i, j) = true;
                end
            end
        catch
            % If fzero fails (e.g., no root found), skip this iteration
            continue;
        end
    end
end

% Plot the regions
figure;
hold on;

% Fill the area that meets the conditions for region1 in gray
contourf(B_n, I_c, region1, [1 1], 'LineColor', 'none', 'FaceColor', [0.5 0.5 0.5]);

% Fill the area that meets the conditions for region2 in blue
contourf(B_n, I_c, region2, [1 1], 'LineColor', 'none', 'FaceColor', 'blue');
grid on;
box on;
set(gca,'LineWidth',2,'FontSize',13);
set(gca, 'GridLineStyle' ,':','linewidth',1');
xlabel('$B_n$','Interpreter','LaTex','FontSize',13);
ylabel('$I_c$','Interpreter','LaTex','FontSize',13);
legend('$E_0^*$ and $E_1^*$','$E_0^*$ and $E^*$','Interpreter','LaTex','FontSize',13)
xlim([0 10]);
ylim([0 10]);
hold off;
