clear
clc

% PARAMETERS
params
p.z_tauchen = 3;
p.A_tauchen = 0.8;

% program
p.Nk = 130; p.maxk = 3; p.mink = 0; maxk = p.maxk; mink = p.mink; % grid capital
p.NA = 5;
p.Nz = 11;
p.NK = 50; 
p.K_range_lo = 0.88; p.K_range_hi = 1.12;
maxK = p.K_range_hi*(maxk-mink); minK = p.K_range_lo*(maxk-mink);
p.alpha_adj = 0.4; % weight on new coefficients
p.T = 2700; % length of simulation


% call function:
use_phi = 0.1;
results = cell(length(use_phi), 1);
for i = 1:length(use_phi)
    p.phi = use_phi(i);
    try
        [coeffs, mean_adjusters, mean_TFPRdisp] = adj_cost_fun(p)
        results{i} = {p, coeffs, mean_adjusters, mean_TFPRdisp};
    catch
        results{i} = 'did not converge';
    end
end
save results results