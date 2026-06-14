function theta_new = thermal_volt(psi, theta_old, theta_a, theta_b, h, mun, Jn_edge, n)
%THERMAL_VOLT Solves the 1D energy balance equation for electron thermal voltage.
%
%   theta_new = thermal_volt(psi, theta_old, theta_a, theta_b, h, mun, Jn_edge, n)
%   computes the updated thermal voltage profile across a 1D spatial grid 
%   using a finite difference scheme with Bernoulli function stabilization.
%
%   Inputs:
%       psi       - Electric potential array evaluated at the grid nodes [V].
%       theta_old - Thermal voltage profile from the previous iteration [V].
%       theta_a   - Dirichlet boundary condition for thermal voltage at x = 0 [V].
%       theta_b   - Dirichlet boundary condition for thermal voltage at x = L [V].
%       h         - Spatial step size (grid spacing) [m].
%       mun       - Electron mobility array evaluated at the nodes [m^2/Vs].
%       Jn_edge   - Electron current density evaluated at the grid edges [A/m^2].
%       n         - Electron density array evaluated at the nodes [m^-3].
%
%   Outputs:
%       theta_new - Updated thermal voltage array including boundaries [V].
%
%   Notes:
%       - This function requires the 'data' script to be available in the 
%         MATLAB path to load global constants (e.g., q, tau_e, Nint, VT).
%       - The linear system is assembled and solved using a sparse tridiagonal matrix.

    data; 
    
    % Edge quantities and transport coefficients
    n_avg = (n(1:end-1) + n(2:end)) / 2;
    theta_avg = (theta_old(1:end-1) + theta_old(2:end)) / 2;
    
    lambda_edge = 5/2 * q * mun .* n_avg .* theta_avg;
    Pe = 5/2 * h * Jn_edge ./ lambda_edge;
    
    [bp, bn] = bern(Pe);
    alpha_i = 3/2 * q * n / tau_e;
    
    % Matrix assembly for the tridiagonal system
    diag_main = lambda_edge(2:end) .* bn(2:end) + lambda_edge(1:end-1) .* bp(1:end-1) + alpha_i(2:end-1) * h^2;
    row_main = (1:Nint)';
    col_main = (1:Nint)';
    
    diag_sub = -lambda_edge(2:end-1) .* bn(2:end-1);
    row_sub = (2:Nint)';
    col_sub = (1:Nint-1)';
    
    diag_super = -lambda_edge(2:end-1) .* bp(2:end-1);
    row_super = (1:Nint-1)';
    col_super = (2:Nint)';
    
    rows = [row_main; row_sub; row_super];
    cols = [col_main; col_sub; col_super];
    vals = [diag_main; diag_sub; diag_super];
    An = sparse(rows, cols, vals, Nint, Nint);
    
    % Source terms and right-hand side vectors
    E_int = -(psi(2:end) - psi(1:end-1)) / h;
    H_i = (Jn_edge(2:end) .* E_int(2:end) + Jn_edge(1:end-1) .* E_int(1:end-1)) / 2;
    
    Fn = h^2 * (H_i + alpha_i(2:end-1) * VT);
    Fn(1)   = Fn(1) + lambda_edge(1) * bn(1) * theta_a;
    Fn(end) = Fn(end) + lambda_edge(end) * bp(end) * theta_b;
    
    % Solve the linear system
    theta_int = An \ Fn;
    theta_new = [theta_a; theta_int; theta_b];
end