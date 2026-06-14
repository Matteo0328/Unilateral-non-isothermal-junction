clear; clc; close all;

%% 1. Initialization and Configuration
% Load physical parameters and mesh
data;

% Simulation flags
PLOTBIAS   = 1; % Plot results per bias step
PLOTTEMP   = 1; % Plot temperature profile
MOBILITY_T = 0; % Temperature-dependent mobility (set only one to 1, both 0 for const)
MOBILITY_E = 1; % Electric field-dependent mobility
PRINT      = 0; % Export PDF plots
ISOTHERMAL = 0; % 1 = Isothermal, 0 = Non-isothermal (Hydrodynamic)

% Initial guess configuration
if ISOTHERMAL == 1
    initial_guess = 1; % Based on previous voltage
elseif ISOTHERMAL == 0
    initial_guess = 2; % Based on isothermal and same voltage case
end

% Solver parameters
if ISOTHERMAL == 1
    tol = 1e-10;
else
    tol = 1e-8;
end
maxit     = 500;
omega_n   = 1;
omega_T   = 1;
omega_psi = 1;

%% 2. Storage Allocation and Data Loading
num_bias = length(Va_list);

% Preallocate solution struct
empty_solution = struct('Va', [], 'psi', [], 'n', [], 'theta', [], ...
    'Jn', [], 'Sn', [], 'Tn', [], 'mun', [], ...
    'JnE', [], 'dSn', []);
solutions = repmat(empty_solution, num_bias, 1);

% Load isothermal baseline if required
if initial_guess == 2
    if isfile('soluzione_iso.mat')
        iso_solutions = load('soluzione_iso.mat');
    else
        error('File "soluzione_iso.mat" not found. Required for initial_guess = 2.');
    end
end

% Create data output directories
cartella_tau = 'Variazione tau';
if ~exist(cartella_tau, 'dir'), mkdir(cartella_tau); end

cartella_Nint = 'Variazione Nint';
if ~exist(cartella_Nint, 'dir'), mkdir(cartella_Nint); end

cartella_mu = 'Variazione mu';
if ~exist(cartella_mu, 'dir'), mkdir(cartella_mu); end

%% 3. Initial guess for the first bias
if initial_guess == 1
    Va = Va_list(1);
    psi_left  = 0;
    psi_right = Va;
    
    % Built-in potential initial guess
    psi = VT*log(ND/ND_high) + Va*x/L;
    psi(1)   = psi_left;
    psi(end) = psi_right;
    
    % Quasi-neutral density initial guess
    n = ND;
    n(1)   = n_left;
    n(end) = n_right;
elseif initial_guess == 2
    Va = iso_solutions.solutions(1).Va;
    psi_left  = 0;
    psi_right = Va;
    
    % Built-in potential initial guess from isothermal
    psi = iso_solutions.solutions(1).psi;
    psi(1)   = psi_left;
    psi(end) = psi_right;
    
    % Quasi-neutral density initial guess from isothermal
    n = iso_solutions.solutions(1).n;
    n(1)   = n_left;
    n(end) = n_right;
end

Tn = TL * ones(Nint+2,1);
theta = VT * ones(Nint+2,1);
mun = mun_0*ones(Nint+1,1); 

%% 4. Bias Continuation Loop
for ibias = 1:num_bias
    Va = Va_list(ibias);
    fprintf('\n=========================================\n');
    fprintf('Solving for Va = %.6g V\n', Va);
    fprintf('=========================================\n');
    
    psi_left  = 0;
    psi_right = Va;
    theta_left = VT;
    theta_right = VT;
    
    if ibias > 1
        if initial_guess == 1
            % Use previous solution as initial guess
            n = solutions(ibias-1).n;
            Va_old = Va_list(ibias-1);
            dVa = Va - Va_old;
            psi = solutions(ibias-1).psi + dVa*x/L;
            psi(1)   = psi_left;
            psi(end) = psi_right;
            theta = solutions(ibias-1).theta;
            mun = solutions(ibias-1).mun;
        elseif initial_guess == 2
            % Use previous isothermal solution as initial guess
            n = iso_solutions.solutions(ibias).n;
            psi = iso_solutions.solutions(ibias).psi;
            psi(1)   = psi_left;
            psi(end) = psi_right;
            theta = VT * ones(Nint+2,1);
            mun = mun_0*ones(Nint+1,1);
        end
    end
    
    err = 1 + tol;
    it  = 0;

    % Gummel Iteration Loop
    while err > tol && it < maxit
        it = it + 1;
        psi_old = psi;
        n_old   = n;
        theta_old = theta;
        
        % A. Continuity equation (Scharfetter-Gummel)
        n_new = continuity_n(psi,theta_old,n_left,n_right,Nint,h,mun);
        n = (1-omega_n)*n_old + omega_n*n_new;
        n(1)   = n_left;
        n(end) = n_right;
        
        % Current at edges
        Jn = electron_current_sg(psi,theta,n,h,mun);
        
        % B. Energy equation (Scharfetter-Gummel)
        if ISOTHERMAL == 0
            theta_new = thermal_volt(psi,theta_old,theta_left,theta_right,h,mun, Jn, n);
            theta = (1-omega_T)*theta_old + omega_T*theta_new;
            theta(1)   = theta_left;
            theta(end) = theta_right;
            Tn = q*theta/kB;
        end
        Sn = energy_flux_sg(n,Jn,theta,mun);
        
        % C. Gummel-linearized Poisson equation
        n_int      = n(2:end-1);
        theta_int  = theta(2:end-1);
        psi_old_in = psi_old(2:end-1);
        ND_int     = ND(2:end-1);
        
        diag_main = 2*epsSi/h^2 + q*n_int./theta_int;
        diag_sub  = -epsSi/h^2 * ones(Nint, 1); 
        diag_sup  = -epsSi/h^2 * ones(Nint, 1);
        
        Mreact = spdiags([diag_sub, diag_main, diag_sup], [-1, 0, 1], Nint, Nint);
        
        b = q*(ND_int - n_int) + q*(n_int./theta_int).*psi_old_in;
        b(1)   = b(1)   + epsSi*psi_left/h^2;
        b(end) = b(end) + epsSi*psi_right/h^2;
        
        psi_new = [psi_left;  Mreact\b; psi_right];
        psi = (1-omega_psi)*psi_old + omega_psi*psi_new;
        psi(1)   = psi_left;
        psi(end) = psi_right;
        
        % D. Mobility update
        mun_new = mun_0 * ones(Nint+1, 1);
        if MOBILITY_T
            Tn_edge = (Tn(1:end-1) + Tn(2:end)) / 2;
            mun_new = mun_new .* (TL ./ Tn_edge).^gamma;
        elseif MOBILITY_E
            E_field = -diff(psi)/h;
            mun_new = mun_new ./ (1 + mun_new .* abs(E_field) / vsat);
        end
        mun = mun_new;
        
        % E. Error computation
        err_psi = norm(psi - psi_old,inf)/max([VT,norm(psi,inf)]);
        err_n   = norm(n - n_old,inf)/max(norm(n,inf),1);
        err_theta   = norm(theta - theta_old,inf)/max(norm(theta,inf),1);
        err = max([err_psi, err_n, err_theta]);
        
        Jrel = (max(Jn)-min(Jn))/max(abs(mean(Jn)),1e-300);
        fprintf('it = %4d | err = %.3e | err_psi = %.3e | err_n = %.3e | err_theta = %.3e | min psi = %.4e | max psi = %.4e \n', ...
            it, err, err_psi, err_n, err_theta, min(psi), max(psi));
            
        if any(theta<=0)
            error("Theta values < 0 detected.");
        end
    end
    
    if it == maxit
        warning('Maximum number of iterations reached at Va = %.4g V',Va);
    end
    
    %% 5. Save Solution
    solutions(ibias).Va     = Va;
    solutions(ibias).psi    = psi;
    solutions(ibias).n      = n;
    solutions(ibias).theta  = theta;
    solutions(ibias).Jn     = electron_current_sg(psi,theta,n,h,mun);
    solutions(ibias).Sn     = energy_flux_sg(n,Jn,theta,mun);
    solutions(ibias).Tn     = Tn;
    solutions(ibias).mun    = mun;
    
    E_int = - (psi(2:end)-psi(1:end-1))/h;
    H_i = (Jn(2:end).*E_int(2:end) + Jn(1:end-1).*E_int(1:end-1))/2;
    H_i = [0; H_i; 0];
    solutions(ibias).JnE = H_i;
    
    dSn = (solutions(ibias).Sn(2:end)-solutions(ibias).Sn(1:end-1))/h;
    solutions(ibias).dSn    = dSn;
    
    %% 6. Plot Current Solution
    if PLOTBIAS  && ISOTHERMAL == 0
        figure(ibias); clf;
        t = tiledlayout(7,1, 'TileSpacing', 'compact', 'Padding', 'compact');
        
        % Electron density
        nexttile;
        semilogy(x*1e6,n,'b',x*1e6, iso_solutions.solutions(ibias).n,'k--','LineWidth',1.4); 
        hold on;
        xline(LC*1e6, '--k', 'LineWidth', 0.7);
        xline((L-LC)*1e6, '--k', 'LineWidth', 0.7);
        grid on;
        xlabel('x [\mum]');
        ylabel('density [m^{-3}]');
        title('Electron density');
        
        % Electrostatic potential
        nexttile;
        p1 = plot(x*1e6,psi,'b', x*1e6, iso_solutions.solutions(ibias).psi,'k--','LineWidth',1.4); 
        hold on;
        xline(LC*1e6, '--k', 'LineWidth', 0.7);
        xline((L-LC)*1e6, '--k', 'LineWidth', 0.7);
        grid on;
        xlabel('x [\mum]');
        ylabel('\psi [V]');
        title('Electrostatic potential');
        
        % SG electron current density
        nexttile;
        Jn_plot = solutions(ibias).Jn;
        plot(x_edge*1e6,round(Jn_plot,2),'b',x_edge*1e6,round(iso_solutions.solutions(ibias).Jn),'k--','LineWidth',1.4);
        hold on;
        xline(LC*1e6, '--k', 'LineWidth', 0.7);
        xline((L-LC)*1e6, '--k', 'LineWidth', 0.7);
        grid on;
        xlabel('x [\mum]');
        ylabel('J_n [A/m^2]');
        title('SG electron current density');
        ytickformat('%.2f');
        
        % Quasi-neutrality ratio
        nexttile;
        plot(x*1e6,n./ND,'b',x*1e6,iso_solutions.solutions(ibias).n./ND,'k--','LineWidth',1.4);
        hold on;
        xline(LC*1e6, '--k', 'LineWidth', 0.7);
        xline((L-LC)*1e6, '--k', 'LineWidth', 0.7);
        grid on;
        xlabel('x [\mum]');
        ylabel('n/N_D');
        title('Quasi-neutrality ratio');
        
        % Space charge density
        nexttile;
        plot(x*1e6,q*(ND-n),'b',x*1e6,q*(ND - iso_solutions.solutions(ibias).n),'k--','LineWidth',1.4);
        hold on;
        xline(LC*1e6, '--k', 'LineWidth', 0.7);
        xline((L-LC)*1e6, '--k', 'LineWidth', 0.7);
        grid on;
        xlabel('x [\mum]');
        ylabel('\rho = q(N_D-n) [C/m^3]');
        title('Space charge density');
        
        % Temperature
        nexttile;
        plot(x*1e6, round(q*theta/kB, 2), 'b', x*1e6,TL*ones(size(x)),'k--', 'LineWidth', 1.4);
        hold on;
        xline(LC*1e6, '--k', 'LineWidth', 0.7);
        xline((L-LC)*1e6, '--k', 'LineWidth', 0.7);
        grid on;
        xlabel('x [\mum]');
        ylabel('T [K]');
        title('Temperature');
        
        % Energy terms
        nexttile;
        p2 = plot(x*1e6,H_i,'-r', x*1e6,3/2*q*n.*(theta-VT)/tau_e,'-c', x(2:end-1)*1e6,round(dSn,2),'-g','LineWidth',1.4);
        for i = 1:length(p2)
            p2(i).Color(4) = 0.7;
        end
        hold on;
        xline(LC*1e6, '--k', 'LineWidth', 0.7);
        xline((L-LC)*1e6, '--k', 'LineWidth', 0.7);
        grid on;
        xlabel('x [\mum]');
        ylabel('[W/m^3]');
        title('Energy terms');
        
        % Global title and legend
        if MOBILITY_E == 0 && MOBILITY_T == 0
            title(t, sprintf('Va = %.3g V, \\tau_E = %.3g ps,N = %.3g, \\mu = const',Va , tau_e*1e12, N),'FontWeight', 'bold', 'FontSize', 16);
        else
            title(t, sprintf('Va = %.3g V, \\tau_E = %.3g ps, N = %.3g',Va , tau_e*1e12, N),'FontWeight', 'bold', 'FontSize', 16);
        end
    
        lgd = legend([p1(1), p1(2), p2(1), p2(2), p2(3)], {'Non-isothermal', ...
            'Isothermal','$J_n\cdot E$', '$\frac{3}{2}qn\frac{\theta -V_T}{\tau_E}$', ...
            '$\frac{dS_n}{dx}$'}, 'Orientation', 'horizontal', 'Interpreter', 'latex', 'FontSize', 14);
        lgd.Layout.Tile = 'south';
        drawnow;
    end
    
    %% 7. Final Residual Diagnostics
    rho_int = q*(ND(2:end-1) - n(2:end-1));
    lhsP = -epsSi*(psi(3:end)-2*psi(2:end-1)+psi(1:end-2))/h^2;
    rhsP = rho_int;
    resP = lhsP - rhsP;
    
    resC_flux = diff(Jn);
    resC_div  = resC_flux/h;
    
    Jscale = q*Dn*ND_high/L;
    Pscale = q*ND_high;
    Cscale = Jscale/L;
    
    meanJn = mean(Jn);
    absSpreadJn = max(Jn) - min(Jn);
    psi_iso = VT*log(ND/ND_high);
    
    % Energy residual
    alpha_i = 3/2 * q * n(2:end-1) /tau_e;
    resE = (Sn(2:end)-Sn(1:end-1))/h - (H_i(2:end-1) - alpha_i.*(theta(2:end-1)-VT));
    resE_norm = norm(resE,inf)/(max(abs(H_i(2:end-1))+abs(alpha_i.*(theta(2:end-1) - VT)) + eps));
    
    % Global energy balance
    lintJnE = trapz(x,H_i);
    rintJnE = Sn(end) - Sn(1) + 3/2*q*trapz(x,n.*(theta-VT)/tau_e);
    
    fprintf('\nFinal diagnostics for Va = %.6e V\n', Va);
    
    % Potential
    fprintf('min(psi)                  = %.6e V\n', min(psi));
    fprintf('max(psi)                  = %.6e V\n', max(psi));
    if ibias == 1
        fprintf('min(psi-psi_iso)          = %.6e V\n', min(psi-psi_iso));
        fprintf('max(psi-psi_iso)          = %.6e V\n', max(psi-psi_iso));
        fprintf('mean(psi-psi_iso)         = %.6e V\n', mean(psi-psi_iso));
        fprintf('Rel. psi error            = %.6e V\n', max(abs(psi-psi_iso))/max(abs(psi_iso)));
    end
    
    % Density
    fprintf('min(n)                    = %.6e m^-3\n', min(n));
    fprintf('max(n)                    = %.6e m^-3\n', max(n));
    fprintf('min(n/ND)                 = %.6e\n', min(n./ND));
    fprintf('max(n/ND)                 = %.6e\n', max(n./ND));
    fprintf('mean(n/ND)                = %.6e\n', mean(n./ND));
    
    % Current
    fprintf('mean(Jn)                  = %.6e A/m^2\n', meanJn);
    fprintf('min(Jn)                   = %.6e A/m^2\n', min(Jn));
    fprintf('max(Jn)                   = %.6e A/m^2\n', max(Jn));
    fprintf('current absolute spread   = %.6e A/m^2\n', absSpreadJn);
    fprintf('spread / Jscale           = %.6e\n', absSpreadJn/Jscale);
    if abs(meanJn) > 1e-12*Jscale
        fprintf('spread / |mean(Jn)|       = %.6e\n', absSpreadJn/abs(meanJn));
    else
        fprintf('spread / |mean(Jn)|       = not meaningful, mean current nearly zero\n');
    end
    
    % Temperature
    fprintf('mean(Tn)                  = %.6e K\n', mean(Tn));
    fprintf('min(Tn)                   = %.6e K\n', min(Tn));
    fprintf('max(Tn)                   = %.6e K\n', max(Tn));
    fprintf('Max rel. temp. increment  = %.6e\n', max((Tn-TL)/TL));
    
    % Poisson residual
    fprintf('Poisson residual inf      = %.6e C/m^3\n', norm(resP,inf));
    fprintf('Poisson residual / Pscale = %.6e\n', norm(resP,inf)/Pscale);
    fprintf('Poisson residual / max(abs(rho_int)) = %.6e\n', norm(resP,inf)/max(norm(rho_int,inf),1));
    
    % Continuity residual
    fprintf('Continuity flux residual  = %.6e A/m^2\n', norm(resC_flux,inf));
    fprintf('Continuity div residual   = %.6e A/m^3\n', norm(resC_div,inf));
    fprintf('Continuity div / Cscale   = %.6e\n', norm(resC_div,inf)/Cscale);
    
    % Energy residual
    fprintf('Energy residual norm        = %.6e \n', resE_norm);
    
    % Density residual
    fprintf('min(rho)                  = %.6e C/m^3\n', min(q*(ND-n)));
    fprintf('max(rho)                  = %.6e C/m^3\n', max(q*(ND-n)));
    fprintf('max(abs(rho))/Pscale      = %.6e\n', max(abs(q*(ND-n)))/Pscale);
    
    % Energy imbalance
    fprintf('Energy imbalance          = %.6e W/m^2\n', abs(lintJnE-rintJnE));
    
    % Drift velocity diagnostics
    vsat = 1e5;
    n_edge = 0.5*(n(1:end-1) + n(2:end));
    u_edge = -Jn./(q*n_edge);
    
    fprintf('max(|u_n|)              = %.6e m/s\n', max(abs(u_edge)));
    fprintf('max(|u_n|)/vsat         = %.6e\n', max(abs(u_edge))/vsat);
    
    Nref = ND_high;
    psiF = psi - VT*log(n/Nref);
    Fqf = -diff(psiF)/h;
    u_qf = mun.*abs(Fqf);
    
    fprintf('max(mu*|grad psiF|)     = %.6e m/s\n', max(u_qf));
    fprintf('max(mu*|grad psiF|)/vsat = %.6e\n', max(u_qf)/vsat);
    
    E = -diff(psi)/h;
    v_drift_only = mun.*abs(E);
    fprintf('max(mu*|E|)/vsat        = %.6e  <-- drift-only, includes built-in field\n', max(v_drift_only)/vsat);
end

%% 8. Global Drift Velocity Diagnostics
for ibias = 1:length(solutions)
    Va_loc  = solutions(ibias).Va;
    psi_loc = solutions(ibias).psi;
    n_loc   = solutions(ibias).n;
    Jn_loc  = solutions(ibias).Jn;
    Tn_loc  = solutions(ibias).Tn;
    mun_loc = solutions(ibias).mun;
    
    Jscale = q*Dn*ND_high/L;
    meanJn = mean(Jn_loc);
    absSpreadJn = max(Jn_loc) - min(Jn_loc);
    
    n_edge = 0.5*(n_loc(1:end-1)+n_loc(2:end));
    u_edge = -Jn_loc./(q*n_edge);
    
    Nref = ND_high;
    psiF = psi_loc - VT*log(n_loc/Nref);
    Fqf = -diff(psiF)/h;
    
    fprintf(['Va = %.3f V | max(|u|)/vsat = %.3e | max(mu|grad psiF|)/vsat = %.3e ' ...
        '| max(T_n) = %.3e | spread(Jn)/Jscale = %.3e | spread(Jn)/|Jmean| = %.3e\n'],  ...
        Va_loc, max(abs(u_edge))/vsat, max(mun_loc.*abs(Fqf))/vsat, max(Tn_loc), absSpreadJn/Jscale, absSpreadJn/abs(meanJn));
end

%% 9. I-V Characteristics Plot
Va_values = arrayfun(@(s) s.Va, solutions);
J_values  = arrayfun(@(s) mean(s.Jn), solutions);

if ISOTHERMAL == 0
    Va_values_iso = arrayfun(@(s) s.Va, iso_solutions.solutions);
    J_values_iso  = arrayfun(@(s) mean(s.Jn), iso_solutions.solutions);
end

figure;
plot(Va_values,abs(J_values),'bo-','LineWidth',1.4);
if ISOTHERMAL == 0
    hold on;
    plot(Va_values_iso,abs(J_values_iso),'ro-','LineWidth',1.4);
end
grid on;
xlabel('V_a [V]');
ylabel('mean J_n [A/m^2]');
title('Current-voltage characteristic');
if ISOTHERMAL == 0
    legend('Non-isothermal', 'Isothermal');
end

%% 10. Temperature Plot
if PLOTTEMP
    figure;
    hold on;
    num_bias = length(solutions);
    h_linee = gobjects(num_bias, 1);
    nomi_legenda = cell(num_bias, 1);
    matrice_colori = turbo(num_bias);
    
    for ibias = 1:num_bias
        colore_temp = matrice_colori(ibias, :);
        h_linee(ibias) = plot(x*1e6,solutions(ibias).Tn, 'Color', colore_temp, 'LineWidth', 1.4); 
        nomi_legenda{ibias} = sprintf('V_a = %.2f V', solutions(ibias).Va);
    end
    grid on;
    xlabel('x [\mum]');
    ylabel('Temperature T [K]');
    title('Temperature');
    legend(h_linee, nomi_legenda);
    drawnow;
end

%% 11. Energy Balance Printout
for ibias = 1 : length(Va_list)
    Va_loc  = solutions(ibias).Va;
    n_loc   = solutions(ibias).n;
    Sn_loc  = solutions(ibias).Sn;
    H_i_loc = solutions(ibias).JnE;
    theta_loc = solutions(ibias).theta;
    
    lintJnE = trapz(x,H_i_loc);
    rintJnE = Sn_loc(end) - Sn_loc(1) + 3/2*q*trapz(x,n_loc.*(theta_loc-VT)/tau_e);
    resJnE = lintJnE - rintJnE;
    
    fprintf('Va = %.3f V | Residuo del bilancio energetico = %.3e \n', Va_loc, resJnE);
end

%% 12. Data and Output Saving
if ISOTHERMAL == 1
    save("soluzione_iso.mat","solutions");
    if N == 1000
        nomeFileDinamico = sprintf('Soluzioni_mu_iso.mat');
        percorsoCompleto = fullfile(cartella_mu, nomeFileDinamico);
        save(percorsoCompleto, "solutions");
    end
elseif ISOTHERMAL == 0
    if N == 1000
        if MOBILITY_T == 1 && MOBILITY_E == 0
            nomeFileDinamico = sprintf('Soluzioni_mu_T_tau=%.2f.mat', tau_e*1e12*100);
            percorsoCompleto = fullfile(cartella_mu, nomeFileDinamico);
            save(percorsoCompleto, "solutions");
        elseif MOBILITY_T == 0 && MOBILITY_E == 1
            nomeFileDinamico = sprintf('Soluzioni_mu_E_tau=%.2f.mat', tau_e*1e12*100);
            percorsoCompleto = fullfile(cartella_mu, nomeFileDinamico);
            save(percorsoCompleto, "solutions");
        elseif MOBILITY_T == 0 && MOBILITY_E == 0
            nomeFileDinamico = sprintf('Soluzioni_mu_0.mat');
            percorsoCompleto = fullfile(cartella_mu, nomeFileDinamico);
            save(percorsoCompleto, "solutions");
            
            nomeFileDinamico = sprintf('Soluzioni_Nint=%d.mat', Nint);
            percorsoCompleto = fullfile(cartella_Nint, nomeFileDinamico);
            save(percorsoCompleto, "solutions");
            
            nomeFileDinamico = sprintf('Soluzioni_tau=%.2f.mat', tau_e*1e12);
            percorsoCompleto = fullfile(cartella_tau, nomeFileDinamico);
            save(percorsoCompleto, "solutions");
        end
    elseif tau_e == 0.2e-12 && N ~= 1000 % update of data in case of a Nint sweep
        nomeFileDinamico = sprintf('Soluzioni_Nint=%d.mat', Nint);
        percorsoCompleto = fullfile(cartella_Nint, nomeFileDinamico);
        save(percorsoCompleto, "solutions");
    end
end

if ISOTHERMAL == 0 && PRINT == 1 && PLOTBIAS == 1 && PLOTTEMP == 1
    exportPDF('Va=0.pdf',figure(1), 'Va=0.05.pdf',figure(2), 'Va=0.1.pdf',figure(3), ...
            'Va=0.2.pdf',figure(4), 'Va=0.3.pdf',figure(5), 'Va=0.4.pdf',figure(6), ...
            'Va=0.5.pdf',figure(7), 'A4Scale', true, 'Margin', 1.5, 'Orientation', ...
            'portrait', 'KeepAspect', false, Folder="V_sweep");
    exportPDF('Temp.pdf',figure(9), Folder="V_sweep" );
end