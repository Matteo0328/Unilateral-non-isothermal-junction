function n = continuity_n(psi, theta, na, nb, Nint, h, mun)
%CONTINUITY_N Solves the 1D steady-state electron continuity equation.
%
%   n = continuity_n(psi, theta, na, nb, Nint, h, mun) computes the 
%   electron density profile across a 1D spatial grid using the 
%   Scharfetter-Gummel finite difference discretization scheme.
%
%   Inputs:
%       psi   - Electric potential array evaluated at the nodes [V].
%       theta - Thermal voltage profile array [V].
%       na    - Dirichlet boundary condition for electron density at x = 0.
%       nb    - Dirichlet boundary condition for electron density at x = L.
%       Nint  - Number of interior grid nodes.
%       h     - Spatial step size (grid spacing) [m].
%       mun   - Electron mobility array evaluated at the nodes [m^2/Vs].
%
%   Outputs:
%       n     - Updated electron density array including boundaries.
%
%   Notes:
%       - This function requires the 'data' script to load necessary globals.
%       - Requires the 'logaritmic_mean' and 'bern' auxiliary functions.

    data; 
    
    % Compute thermal voltage and potential gradients at grid edges
    theta_med = logaritmic_mean(theta(2:end), theta(1:end-1));
    eta = (diff(psi) - diff(theta)) ./ theta_med;
    
    % Evaluate Bernoulli functions for Scharfetter-Gummel stabilization
    [bp, bn] = bern(eta);
    
    % Assemble the tridiagonal matrix for the interior nodes
    diag_main = -(mun(2:end) .* theta_med(2:end) .* bn(2:end) + mun(1:end-1) .* theta_med(1:end-1) .* bp(1:end-1));
    row_main = (1:Nint)';
    col_main = (1:Nint)';
    
    diag_sub = mun(2:end-1) .* theta_med(2:end-1) .* bn(2:end-1);
    row_sub = (2:Nint)';
    col_sub = (1:Nint-1)';
    
    diag_super = mun(2:end-1) .* theta_med(2:end-1) .* bp(2:end-1);
    row_super = (1:Nint-1)';
    col_super = (2:Nint)';
    
    rows = [row_main; row_sub; row_super];
    cols = [col_main; col_sub; col_super];
    vals = [diag_main; diag_sub; diag_super];
    An = sparse(rows, cols, vals, Nint, Nint);
    
    % Apply Dirichlet boundary conditions to the load vector
    Fn = sparse(Nint, 1);
    Fn(1)   = -mun(1)   * theta_med(1)   * bn(1)   * na;
    Fn(end) = -mun(end) * theta_med(end) * bp(end) * nb;
    
    % Solve the linear system and append boundary nodes
    n_int = An \ Fn;
    n = [na; n_int; nb];
    
    % (Optional/Commented: Edge current density calculation)
    % Jn_edge = (q*mun)/h * theta_med .* (n(2:end).*bp - n(1:end-1).*bn);
end