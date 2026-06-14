% Initialize workspace and load spatial grid parameters
clear; clc; close all;
addpath('..');
data; 

% Load simulation results for various energy relaxation times (tau_E)
s005 = load("Soluzioni_tau=0.05.mat");
s01  = load("Soluzioni_tau=0.10.mat");
s02  = load("Soluzioni_tau=0.20.mat");
s05  = load("Soluzioni_tau=0.50.mat");
s1   = load("Soluzioni_tau=1.00.mat");

% Loop through each applied bias voltage
for ibias = 1 : length(Va_list)
    figure(ibias); clf;
    
    % Extract electron density (n)
    n005 = s005.solutions(ibias).n;
    n01  = s01.solutions(ibias).n;
    n02  = s02.solutions(ibias).n;
    n05  = s05.solutions(ibias).n;
    n1   = s1.solutions(ibias).n;
    
    % Extract energy flux (dSn)
    dSn005 = s005.solutions(ibias).dSn;
    dSn01  = s01.solutions(ibias).dSn;
    dSn02  = s02.solutions(ibias).dSn;
    dSn05  = s05.solutions(ibias).dSn;
    dSn1   = s1.solutions(ibias).dSn;
    
    % Extract electron temperature (Tn)
    Tn005 = s005.solutions(ibias).Tn;
    Tn01  = s01.solutions(ibias).Tn;
    Tn02  = s02.solutions(ibias).Tn;
    Tn05  = s05.solutions(ibias).Tn;
    Tn1   = s1.solutions(ibias).Tn;
    
    % Extract current density (Jn)
    Jn005 = s005.solutions(ibias).Jn;
    Jn01  = s01.solutions(ibias).Jn;
    Jn02  = s02.solutions(ibias).Jn;
    Jn05  = s05.solutions(ibias).Jn;
    Jn1   = s1.solutions(ibias).Jn;
    
    % Extract thermal voltage equivalents (theta)
    theta005 = s005.solutions(ibias).theta;
    theta01  = s01.solutions(ibias).theta;
    theta02  = s02.solutions(ibias).theta;
    theta05  = s05.solutions(ibias).theta;
    theta1   = s1.solutions(ibias).theta;
    
    % --- Plot 1: Electron Temperature ---
    subplot(3,1,1)
    plot(x*1e6, round(Tn005,2), x*1e6, round(Tn01,2), x*1e6, round(Tn02,2), ...
         x*1e6, round(Tn05,2), x*1e6, round(Tn1,2), 'LineWidth', 1.4)
    legend('\tau_E = 0.05', '\tau_E = 0.1', '\tau_E = 0.2', ...
           '\tau_E = 0.5', '\tau_E = 1');
    grid on;
    xlabel('x [\mum]');
    ylabel('Temperature [K]');
    title('Temperature');
    
    % Calculate Lattice Energy Relaxation Term [W/m^3]
    energy005 = 3/2 * q * n005 .* (theta005 - VT) / (0.05 * 1e-12);
    energy01  = 3/2 * q * n01  .* (theta01  - VT) / (0.1  * 1e-12);
    energy02  = 3/2 * q * n02  .* (theta02  - VT) / (0.2  * 1e-12);
    energy05  = 3/2 * q * n05  .* (theta05  - VT) / (0.5  * 1e-12);
    energy1   = 3/2 * q * n1   .* (theta1   - VT) / (1.0  * 1e-12);
    
    % --- Plot 2: Energy Relaxation Rate ---
    subplot(3,1,2)
    plot(x*1e6, energy005, x*1e6, energy01, x*1e6, energy02, x*1e6, energy05, ...
         x*1e6, energy1, 'LineWidth', 1.4)
    grid on;
    xlabel('x [\mum]');
    ylabel('Energy Relaxation Rate [W/m^3]'); 
    title('Lattice Energy Relaxation Term');
    legend('\tau_E = 0.05 ps', '\tau_E = 0.1 ps', '\tau_E = 0.2 ps', ...
           '\tau_E = 0.5 ps', '\tau_E = 1.0 ps');
           
    % --- Plot 3: Current Density ---
    subplot(3,1,3)
    plot(x_edge*1e6, round(Jn005,2), x_edge*1e6, round(Jn01,2), x_edge*1e6, round(Jn02,2), ...
         x_edge*1e6, round(Jn05,2), x_edge*1e6, round(Jn1,2), 'LineWidth', 1.4)
    legend('\tau_E = 0.05', '\tau_E = 0.1', '\tau_E = 0.2', ...
           '\tau_E = 0.5', '\tau_E = 1');
    grid on;
    xlabel('x [\mum]');
    ylabel('Current Density [A/m^2]');
    title('Current Density');
    
    % Print peak values to the command window
    fprintf("Va = %.3g\n", Va_list(ibias))
    fprintf("     |     0.05    |     0.1     |      0.2     |     0.5     |      1  \n");
    fprintf("  Tn |  %.3e  |  %.3e  |  %.3e   |  %.3e  |  %.3e  \n", ...
        max(Tn005), max(Tn01), max(Tn02), max(Tn05), max(Tn1));
    fprintf("  Jn |  %.3e |  %.3e |  %.3e  |  %.3e |  %.3e  \n", ...
        max(Jn005), max(Jn01), max(Jn02), max(Jn05), max(Jn1));
    fprintf(" Energy |  %.3e |  %.3e |  %.3e  |  %.3e |  %.3e  \n", ...
        max(energy005), max(energy01), max(energy02), max(energy05), max(energy1));
end

% Export figures to PDF (Currently commented out)
%exportPDF('Va005.pdf',figure(2),'Va04.pdf',figure(6), 'A4Scale', true, 'Orientation', ...
%        'portrait', 'KeepAspect', false);