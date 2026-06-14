clear; clc; close all;
addpath('..');
data; % Load physical constants (L, Va_list, etc.)

%% 1. Dynamic loading of all Nint files
% Find all Nint simulation files in the current directory
file_list = dir('Soluzioni_Nint=*.mat');
if isempty(file_list)
    error('No files found matching "Soluzioni_Nint=*.mat". Please check your directory.');
end

num_files = length(file_list);
Nint_list = zeros(num_files, 1);
solutions_all = cell(num_files, 1);

% Extract Nint values from filenames and load data
for k = 1:num_files
    filename = file_list(k).name;
    tokens = regexp(filename, 'Nint=(\d+)', 'tokens');
    Nint_list(k) = str2double(tokens{1}{1});
    
    tmp = load(filename);
    solutions_all{k} = tmp.solutions;
end

% Sort loaded files in ascending order based on Nint
[Nint_list, sort_idx] = sort(Nint_list);
solutions_all = solutions_all(sort_idx);

fprintf('Files found and sorted by Nint: %s\n\n', num2str(Nint_list'));


%% 2. Profile plotting (Restricted to Nint = 400, 800, 1600)
% Locate the indices for the specific grids we want to plot
idx_400  = find(Nint_list == 400);
idx_800  = find(Nint_list == 800);
idx_1600 = find(Nint_list == 1600);

if ~isempty(idx_400) && ~isempty(idx_800) && ~isempty(idx_1600)
    % Reconstruct spatial vectors
    x400  = linspace(0,L,402)'; 
    x800  = linspace(0,L,802)'; 
    x1600 = linspace(0,L,1602)'; 
    x_edge400  = 0.5*(x400(1:end-1)  + x400(2:end));
    x_edge800  = 0.5*(x800(1:end-1)  + x800(2:end));
    x_edge1600 = 0.5*(x1600(1:end-1) + x1600(2:end));
    
    s400_sol  = solutions_all{idx_400};
    s800_sol  = solutions_all{idx_800};
    s1600_sol = solutions_all{idx_1600};

    for ibias = 1 : length(Va_list)
        figure(ibias); clf;
        
        Tn400  = s400_sol(ibias).Tn;
        Tn800  = s800_sol(ibias).Tn;
        Tn1600 = s1600_sol(ibias).Tn;
        
        Jn400  = s400_sol(ibias).Jn;
        Jn800  = s800_sol(ibias).Jn;
        Jn1600 = s1600_sol(ibias).Jn;
        
        % Temperature Plot
        subplot(2,1,1)
        plot(x400*1e6, round(Tn400,5), x800*1e6, round(Tn800,5), x1600*1e6, round(Tn1600,5), 'LineWidth', 1.4);
        legend('N_{int} = 400','N_{int} = 800','N_{int} = 1600');
        grid on;
        xlabel('x [\mum]');
        ylabel('Temperature [K]');
        title('Temperature');
        
        % Current Density Plot
        subplot(2,1,2)
        plot(x_edge400*1e6, round(Jn400,2), x_edge800*1e6, round(Jn800,2), x_edge1600*1e6, round(Jn1600,2), 'LineWidth', 1.4);
        legend('N_{int} = 400','N_{int} = 800','N_{int} = 1600');
        grid on;
        xlabel('x [\mum]');
        ylabel('Current Density [A/m^2]');
        title('Current Density');
        
        fprintf("Va = %.3g V\n", Va_list(ibias));
        fprintf("     |     400     |     800     |      1600 \n");
        fprintf("  Tn |  %.3f    |  %.3f    |  %.3f\n", max(Tn400), max(Tn800), max(Tn1600));
        fprintf("  Jn |  %.3e |  %.3e |  %.3e\n\n", max(Jn400), max(Jn800), max(Jn1600));
    end
else
    warning('Simulations for Nint = 400, 800, and 1600 are not all present. Skipping profile plots.');
end

%% 3. Convergence Analysis (Global Fit Calculation)
% Use the finest grid (last in the sorted list) as the exact reference solution
N_ref = Nint_list(end);
sol_ref = solutions_all{end};
x_ref = linspace(0, L, N_ref+2)';

% Test all other grids against the reference
N_test = Nint_list(1:end-1);
num_test = length(N_test);

if num_test < 1
    error('At least one test grid is required alongside the reference grid to compute the error.');
end

h_vec = L ./ (N_test + 1); % Spatial step sizes
num_bias = length(Va_list);

% Arrays to store the Global Fit convergence order
P_N_GLOBAL     = zeros(num_bias, 1);
P_PSI_GLOBAL   = zeros(num_bias, 1);
P_THETA_GLOBAL = zeros(num_bias, 1);

% Compute errors and convergence rates
for ibias = 1 : num_bias
    err_psi_vec   = zeros(num_test, 1);
    err_n_vec     = zeros(num_test, 1);
    err_theta_vec = zeros(num_test, 1);
    
    % Reference variables
    psi_ref   = sol_ref(ibias).psi;
    n_ref     = sol_ref(ibias).n;
    theta_ref = sol_ref(ibias).theta;
    
    for k = 1:num_test
        N_curr = N_test(k);
        sol_curr = solutions_all{k};
        x_curr = linspace(0, L, N_curr+2)';
        
        % Interpolate reference solution onto the coarser test grid
        psi_ref_interp   = interp1(x_ref, psi_ref, x_curr, 'spline');
        n_ref_interp     = interp1(x_ref, n_ref, x_curr, 'spline');
        theta_ref_interp = interp1(x_ref, theta_ref, x_curr, 'spline');
        
        % Compute maximum spatial error (Infinity Norm)
        err_psi_vec(k)   = norm(sol_curr(ibias).psi - psi_ref_interp, inf);
        err_n_vec(k)     = norm(sol_curr(ibias).n - n_ref_interp, inf);
        err_theta_vec(k) = norm(sol_curr(ibias).theta - theta_ref_interp, inf);
    end
    
    if num_test >= 2
        % Global convergence (Linear regression on log-log scale)
        % polyfit(X, Y, 1) finds the slope (order p) of the best fit line
        fit_n     = polyfit(log10(h_vec), log10(err_n_vec), 1);
        fit_psi   = polyfit(log10(h_vec), log10(err_psi_vec), 1);
        fit_theta = polyfit(log10(h_vec), log10(err_theta_vec), 1);
        
        P_N_GLOBAL(ibias)     = fit_n(1);
        P_PSI_GLOBAL(ibias)   = fit_psi(1);
        P_THETA_GLOBAL(ibias) = fit_theta(1);
    end
end


%% 4. Print Global Convergence Table
fprintf("==============================================================\n");
fprintf(" GLOBAL CONVERGENCE ORDER (Exact reference: Nint = %d)\n", N_ref);
fprintf("==============================================================\n");

if num_test >= 2
    fprintf(" Va (Bias)  |    p_n    |   p_psi   |  p_theta  \n");
    fprintf("--------------------------------------------------------------\n");
    for ibias = 1:num_bias
        fprintf(" %-8.3g V |  %7.5f  |  %7.5f  |  %7.5f  \n", ...
            Va_list(ibias), P_N_GLOBAL(ibias), P_PSI_GLOBAL(ibias), P_THETA_GLOBAL(ibias));
    end
    fprintf("==============================================================\n\n");
else
    fprintf("N/A (Need >= 2 test grids to calculate the global slope)\n\n");
end       