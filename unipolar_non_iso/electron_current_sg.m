function Jn = electron_current_sg(psi, theta, n, h, mun)
%ELECTRON_CURRENT_SG Computes the electron current density at grid edges.
%
%   Jn = electron_current_sg(psi, theta, n, h, mun) calculates the 
%   electron current density across a 1D spatial grid using the 
%   Scharfetter-Gummel finite difference discretization scheme.
%
%   Inputs:
%       psi   - Electric potential array evaluated at the nodes [V].
%       theta - Thermal voltage profile evaluated at the nodes [V].
%       n     - Electron density array evaluated at the nodes [m^-3].
%       h     - Spatial step size (grid spacing) [m].
%       mun   - Electron mobility array evaluated at the edges [m^2/Vs].
%
%   Outputs:
%       Jn    - Electron current density evaluated at the grid edges [A/m^2].
%
%   Notes:
%       - Requires the 'data' script in the path for global constants (e.g., q).
%       - Requires the 'logaritmic_mean' and 'bern' auxiliary functions.

    data; 
    
    % Compute average thermal voltage and driving force parameter (eta) at edges
    theta_med = logaritmic_mean(theta(2:end), theta(1:end-1));
    eta = (diff(psi) - diff(theta)) ./ theta_med;
    
    % Evaluate Bernoulli functions for Scharfetter-Gummel stabilization
    [bp, bn] = bern(eta);
    
    % Calculate the stabilized electron current density
    Jn = q / h * mun .* theta_med .* (n(2:end) .* bp - n(1:end-1) .* bn);
    
end