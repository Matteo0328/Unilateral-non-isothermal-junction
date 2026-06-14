% Basic parameters
% Nint = 800, tau = 0.2ps, gamma = 2.5
clear; clc; close all;
addpath('..');
data; % Load spatial grid variables

% Load simulation results for different mobility models
s_mu_T   = load("Soluzioni_mu_T_tau=020.mat");
s_mu_E   = load("Soluzioni_mu_E_tau=020.mat");
s_mu_0   = load("Soluzioni_mu_0.mat");
s_mu_iso = load("Soluzioni_mu_iso.mat");

% --- Temperature Plots ---
figure; clf;

% Low voltage case (Va roughly 0.1V)
T_T_01 = s_mu_T.solutions(3).Tn;
T_E_01 = s_mu_E.solutions(3).Tn;
T_0_01 = s_mu_0.solutions(3).Tn;
Va_01  = s_mu_T.solutions(3).Va;

subplot(1,2,1);
plot(x*1e6, T_T_01, '-b', x*1e6, T_E_01, '-g', x*1e6, T_0_01, '--k');
grid on;
xlabel('x [\mum]');
ylabel('Temperature [K]');
legend('\mu(T)', '\mu(E)', '\mu = const');
title(sprintf('Va = %.3g', Va_01));

% High voltage case (Va roughly 0.5V)
T_T_05 = s_mu_T.solutions(7).Tn;
T_E_05 = s_mu_E.solutions(7).Tn;
T_0_05 = s_mu_0.solutions(7).Tn;
Va_05  = s_mu_T.solutions(7).Va;

subplot(1,2,2);
plot(x*1e6, T_T_05, '-b', x*1e6, T_E_05, '-g', x*1e6, T_0_05, '--k');
grid on;
xlabel('x [\mum]');
ylabel('Temperature [K]');
legend('\mu(T)', '\mu(E)', '\mu = const');
title(sprintf('Va = %.3g', Va_05));

% --- I-V Curve Plots ---
% Get the applied voltages and average currents for each model
Va_values_mu_T = arrayfun(@(s) s.Va, s_mu_T.solutions);
J_values_mu_T  = arrayfun(@(s) mean(s.Jn), s_mu_T.solutions);

Va_values_mu_E = arrayfun(@(s) s.Va, s_mu_E.solutions);
J_values_mu_E  = arrayfun(@(s) mean(s.Jn), s_mu_E.solutions);

Va_values_mu_0 = arrayfun(@(s) s.Va, s_mu_0.solutions);
J_values_mu_0  = arrayfun(@(s) mean(s.Jn), s_mu_0.solutions);

Va_values_iso  = arrayfun(@(s) s.Va, s_mu_iso.solutions);
J_values_iso   = arrayfun(@(s) mean(s.Jn), s_mu_iso.solutions);

figure;
plot(Va_values_mu_T, abs(J_values_mu_T), 'bo-', ...
     Va_values_mu_E, abs(J_values_mu_E), 'go-', ...
     Va_values_mu_0, abs(J_values_mu_0), 'ko-', ...
     Va_values_iso,  abs(J_values_iso),  'ro-', 'LineWidth', 1.4);
grid on;
xlabel('V_a [V]');
ylabel('mean J_n [A/m^2]');
title('Current-voltage characteristic');
legend('\mu(T)', '\mu(E)', '\mu = const', 'Isothermal');

% --- Mobility Plots ---
% Reload data and plot spatial mobility for the 6th solution
sT02 = load("Soluzioni_mu_T_tau=020.mat");
sE02 = load("Soluzioni_mu_E_tau=020.mat");
s02  = load("Soluzioni_mu_0.mat");

muT02 = sT02.solutions(6).mun;
muE02 = sE02.solutions(6).mun;
mu02  = s02.solutions(6).mun;

figure(3);
plot(x_edge*1e6, muT02, x_edge*1e6, muE02, x_edge*1e6, mu02, 'LineWidth', 1.4);
legend('\mu(T)', '\mu(E)', '\mu_0');
grid on;
xlabel('x [\mum]');
ylabel('Mobility [m^2/Vs]');
title('Mobility');

% Save the mobility plot to a PDF
exportPDF('mobility.pdf', figure(3));