% Basic initialization
clear; clc; close all;
addpath('..');
data; % Load required constants/grid

% Load baseline solution (constant mobility)
tmp_0 = load('Soluzioni_mu_0.mat');
f_0 = fieldnames(tmp_0); 
s_0 = tmp_0.(f_0{1});

% Handle data structure and extract I-V for baseline
if isfield(s_0, 'solutions')
    dati_0 = s_0.solutions;
else
    dati_0 = s_0; 
end
Va_0 = arrayfun(@(s) s.Va, dati_0);
J_0  = arrayfun(@(s) mean(s.Jn), dati_0);

% Find parametric files for mu(T) and mu(E)
soluzioni_mu_T = dir('Soluzioni_mu_T_tau=*.mat');
soluzioni_mu_E = dir('Soluzioni_mu_E_tau=*.mat');

% Sort mu(T) files numerically by tau value
nomi_T = {soluzioni_mu_T.name};
numeri_testo_T = regexp(nomi_T, '(?<=tau=)\d+', 'match');
valori_T = cellfun(@(x) str2double(x{1}) / 100, numeri_testo_T);
[~, indici_T] = sort(valori_T);
soluzioni_mu_T = soluzioni_mu_T(indici_T);

% Sort mu(E) files numerically by tau value
nomi_E = {soluzioni_mu_E.name};
numeri_testo_E = regexp(nomi_E, '(?<=tau=)\d+', 'match');
valori_E = cellfun(@(x) str2double(x{1}) / 100, numeri_testo_E);
[~, indici_E] = sort(valori_E);
soluzioni_mu_E = soluzioni_mu_E(indici_E);

% Setup Figure 1 for mu(E)
fig_E = figure(1);
hold on; grid on;
xlabel('V_a [V]'); ylabel('mean J_n [A/m^2]');
title('Current-voltage characteristic - \mu(E)');
plot(Va_0, abs(J_0), 'ko-', 'LineWidth', 1.5, 'MarkerFaceColor', 'k'); % Baseline plot
leg_E = {'\mu = const'};

% Setup Figure 2 for mu(T)
fig_T = figure(2);
hold on; grid on;
xlabel('V_a [V]'); ylabel('mean J_n [A/m^2]');
title('Current-voltage characteristic - \mu(T)');
plot(Va_0, abs(J_0), 'ko-', 'LineWidth', 1.5, 'MarkerFaceColor', 'k'); % Baseline plot
leg_T = {'\mu = const'};

% Generate distinct colors for the loop
colori = lines(length(soluzioni_mu_E)); 

% Loop through sorted files and plot data
for i = 1:length(soluzioni_mu_E)
    tau_corrente = valori_E(i);
    
    % Load data for current tau
    tmp_E = load(fullfile(soluzioni_mu_E(i).folder, soluzioni_mu_E(i).name));
    f_E = fieldnames(tmp_E); 
    s_E = tmp_E.(f_E{1});
    
    tmp_T = load(fullfile(soluzioni_mu_T(i).folder, soluzioni_mu_T(i).name));
    f_T = fieldnames(tmp_T); 
    s_T = tmp_T.(f_T{1});
    
    % Handle data structures
    if isfield(s_E, 'solutions')
        dati_E = s_E.solutions;
    else
        dati_E = s_E; 
    end
    
    if isfield(s_T, 'solutions')
        dati_T = s_T.solutions;
    else
        dati_T = s_T; 
    end
    
    % Extract voltage and current arrays
    Va_E = arrayfun(@(s) s.Va, dati_E);
    J_E  = arrayfun(@(s) mean(s.Jn), dati_E);
    
    Va_T = arrayfun(@(s) s.Va, dati_T);
    J_T  = arrayfun(@(s) mean(s.Jn), dati_T);
    
    % Add to mu(E) figure
    figure(fig_E); 
    plot(Va_E, abs(J_E), 'o-', 'Color', colori(i,:), 'MarkerFaceColor', colori(i,:), 'LineWidth', 1.4);
    leg_E{end+1} = sprintf('\\tau = %g ps', tau_corrente); 
    
    % Add to mu(T) figure
    figure(fig_T); 
    plot(Va_T, abs(J_T), 'x--', 'Color', colori(i,:), 'LineWidth', 1.4);
    leg_T{end+1} = sprintf('\\tau = %g ps', tau_corrente);
end

% Apply legends and finalize figures
figure(fig_E);
legend(leg_E, 'Location', 'best');
hold off;

figure(fig_T);
legend(leg_T, 'Location', 'best');
hold off;

% Export to PDFs
exportPDF('I-V_mu_E.pdf', figure(fig_E), 'I-V_mu_T.pdf', figure(fig_T));