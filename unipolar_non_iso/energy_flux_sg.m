function Sn = energy_flux_sg(n, Jn_edge, theta, mun)
%ENERGY_FLUX_SG Computes the electron energy flux at the grid edges.
%
%   Sn = energy_flux_sg(n, Jn_edge, theta, mun) calculates the energy 
%   flux across a 1D spatial grid using a Scharfetter-Gummel style 
%   finite difference discretization scheme to ensure numerical stability.
%
%   Inputs:
%       n       - Electron density array evaluated at the nodes [m^-3].
%       Jn_edge - Electron current density evaluated at the edges [A/m^2].
%       theta   - Thermal voltage profile evaluated at the nodes [V].
%       mun     - Electron mobility array evaluated at the nodes [m^2/Vs].
%
%   Outputs:
%       Sn      - Electron energy flux evaluated at the grid edges [W/m^2].
%
%   Notes:
%       - Requires the 'data' script in the path for global constants (e.g., q, h).
%       - Requires the 'bern' function to evaluate the Bernoulli numbers.

    data; 
    
    % Compute average density and thermal voltage at the grid edges
    n_avg = (n(1:end-1) + n(2:end)) / 2;
    theta_avg = (theta(1:end-1) + theta(2:end)) / 2;
    
    % Calculate edge thermal conductivity (lambda) and Péclet number (Pe)
    lambda_edge = 5/2 * q * mun .* n_avg .* theta_avg;
    Pe = 5/2 * h * Jn_edge ./ lambda_edge;
    
    % Evaluate Bernoulli functions for numerical stabilization
    [bp, bn] = bern(Pe);
    
    % Compute the stabilized energy flux
    Sn = lambda_edge .* (theta(1:end-1) .* bn - theta(2:end) .* bp) / h;
    
end