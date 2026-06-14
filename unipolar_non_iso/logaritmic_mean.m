function logm = logaritmic_mean(a, b)
%LOGARITMIC_MEAN Computes the logarithmic mean of two arrays.
%
%   logm = logaritmic_mean(a, b) calculates the logarithmic mean 
%   between corresponding elements of arrays 'a' and 'b'. 
%   It uses a Taylor series expansion by default to prevent numerical 
%   instabilities when 'a' and 'b' are very close in value.
%
%   Inputs:
%       a - First input array.
%       b - Second input array (must be the same size as 'a').
%
%   Outputs:
%       logm - Array containing the logarithmic mean values.

    % Initialize output array
    logm = zeros(size(a));
    taylor = 1; % Toggle for Taylor expansion method
    
    if taylor == 1
        % Variables for Taylor expansion
        alpha = (b + a) ./ 2;
        beta = (b - a) ./ 2;
        gamma = beta ./ alpha;
        
        % Compute using Taylor series approximation
        logm = alpha .* (1 - 1/3 * gamma.^2 - 4/45 * gamma.^4 - 44/945 * gamma.^6); 
    else
        % Identify elements with very small relative differences
        rel_diff = abs(b - a) ./ (a + b);
        small = rel_diff < 1e-8;
        large = ~small;
        
        % Handle small differences to avoid division by zero
        if any(small)
            logm(small) = a(small); 
            % Alternatively, use arithmetic mean: 
            % logm(small) = (a(small) + b(small)) ./ 2; 
        end
        
        % Compute standard logarithmic mean for the rest
        if any(large)
            logm(large) = (b(large) - a(large)) ./ (log(b(large)) - log(a(large)));
        end
    end
end