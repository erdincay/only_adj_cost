function Yhats = guess_from_coeffs(coeffs, K, p)
%guess_from_coeffs predict aggregate variables from regression coefficients

% make regressors - one for each aggregate state
new_regs = zeros(p.NK * p.NA, p.NA + 1);
for iA = 1:p.NA
    new_regs((iA - 1) * p.NK + 1:iA * p.NK, iA) = 1; % dummies for tech state
end
new_regs(:,p.NA + 1) = repmat(log(K)', [p.NA 1]); % values for aggregate K

% 'prediction' of aggregate variables
Yhats = new_regs * coeffs;