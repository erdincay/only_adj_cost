clear
clc

tic;
% PARAMETERS
p.beta = 0.98; % discount factor
p.delta = 0.16; % depreciation
p.alpha = 0.33; % production function: y = k^alpha*(zAl)^(1-alpha)
p.rho_A = 0.9; % persistence in A
p.sig_A = 0.03; % stddev in shocks to A
p.gA = 0.01; % growth rate in A (use to ensure eventual adjustment)
p.rho_z = 0.7; % persistence in z
p.sig_z = 0.4; % stddev in shocks to z
p.sigma = 2; % elasticity of substitution between goods
p.phi = 0.01; % fixed adjustment cost (0 = frictionless)
p.psi = 1.05; % CRRA coefficient on household utility

p.sig_Unc = 0.3;
p.sig_A = 0;

p.sigexp = (p.sigma-1)/p.sigma;
p.z_tauchen = 2;
p.A_tauchen = 0.8;

% program
p.Nk = 80; p.maxk = 1.3; p.mink = 0; maxk = p.maxk; mink = p.mink; % grid capital
p.NA = 5;
p.Nz = 7;
p.NK = 50; 
p.K_range_lo = 0.92; p.K_range_hi = 1.12;
maxK = p.K_range_hi*(maxk-mink); minK = p.K_range_lo*(maxk-mink);
alpha_adj_init = 0.3;
p.alpha_adj = alpha_adj_init; % weight on new coefficients
p.T = 1700; % length of simulation

% flags
flags.new_calc = 1;
flags.new_shocks = 0;
flags.matrices_to_update_Mu = 0;

% Set up grids
k = mink:(maxk-mink)/(p.Nk):maxk; k(1)=[];
K = minK:(maxK-minK)/(p.NK):maxK; K(1)=[];
[Unc, pi.Unc] = TAUCHEN(p.NA, p.rho_A, p.sig_Unc, p.A_tauchen); Unc = exp(-Unc');
[A, pi.A] = TAUCHEN(p.NA, p.rho_A, p.sig_A, p.A_tauchen); A = exp(A');
if p.sig_A == 0; pi.A = pi.Unc; end
z = zeros(p.NA, p.Nz);
pi.z = zeros(p.Nz, p.Nz, p.NA);
for iA = 1:p.NA
    [z(iA,:), pi.z(:,:,iA)] = TAUCHEN(p.Nz, p.rho_z, p.sig_z * Unc(iA), p.z_tauchen / Unc(iA));
end
z = exp(z(1,:));
[pi.z, ~, ~] = TAUCHEN_Bloom(p.Nz, z(1), z(end), p.NA, p.sig_z * Unc, p.rho_z, 2, 0, 1);

% construct vector of aggregate shocks:
Astate = aggregate_shocks(pi, p, flags);

% get guesses from economy without aggregate shocks
noA = adj_cost_noA(p);
display('Initial values fetched');

% Set up grids not relying on k
% V = zeros(p.Nz, p.Nk, p.NA, p.NK);
V = repmat(noA.V, [1 1 p.NA p.NK]);
maxk = noA.maxk;

zgrid = repmat(z', [1 p.Nk p.NA p.NK]);
kgrid = repmat(k, [p.Nz 1 p.NA p.NK]);
Agrid = repmat(permute(A, [1 3 2]), [p.Nz p.Nk 1 p.NK]);

% Initial values for coefficients
coeffs = zeros(p.NA + 1, 4); N_coeff = 5; % K, C, P, I, w
coefftol_hist = zeros(1,10);
coeffs(1:p.NA,1) = noA.K_agg; % agg_K;
coeffs(1:p.NA,2) = noA.C_agg; % C;
coeffs(1:p.NA,3) = noA.I_agg; % I;
coeffs(1:p.NA,4) = noA.wguess; % w;

% R_old = zeros(p.NA,p.NK) + 1 / p.beta;
% R_flag = 1;
% alpha_R = 0.1;

% Start coefficients loop
coefftol = 1; coeffct = 0;
while coefftol > 0.004 && coeffct < 90

% Set up grids relying on k
k = mink:(maxk-mink)/(p.Nk):maxk; k(1)=[];
maxK = p.K_range_hi * noA.K_agg; minK = p.K_range_lo * noA.K_agg;
K = minK:(maxK-minK)/(p.NK):maxK; K(1)=[];

kgrid = repmat(k, [p.Nz 1 p.NA p.NK]);

% get guesses from coefficients
% guess = guess_agg_vars_from_coeffs(p, coeffs, K, pi);
Y_hat = guess_from_coeffs(coeffs, K, p);
guess.Kprime = reshape(Y_hat(:,1), [p.NK, p.NA])';
guess.C = reshape(Y_hat(:,2), [p.NK, p.NA])';
guess.I = reshape(Y_hat(:,3), [p.NK, p.NA])';
guess.w = reshape(Y_hat(:,4), [p.NK, p.NA])';

% R = zeros(p.NA, p.NK);
mu_C = guess.C .^ (-p.psi); % marginal utility of consumption
if coeffct == 0
    P_mu = mu_C;
elseif coefftol < 0.2
    P_mu = 0.05 * mu_C + (1 - 0.05) * P_mu;
    P_mu_tol = max(abs(P_mu(:) - mu_C(:)));
    display(P_mu_tol)
end

Igrid = repmat(permute(guess.I, [3 4 1 2]), [p.Nz p.Nk 1 1]);
Pgrid = ones(size(Igrid));
wgrid = repmat(permute(guess.w, [3 4 1 2]), [p.Nz p.Nk 1 1]);
% Rgrid = repmat(permute(R, [3 4 1 2]), [p.Nz p.Nk 1 1]);
P_mugrid = repmat(permute(P_mu, [3 4 1 2]), [p.Nz p.Nk 1 1]);

% calculate in-period returns for all states
denom = p.alpha*(p.sigma-1)+1;
labor = ( ((p.sigma-1)/p.sigma*(1-p.alpha)./wgrid).^(p.sigma) .* (Igrid.*Pgrid.^(p.sigma-1)) .* ...
    kgrid.^(p.alpha*(p.sigma-1)) .* ...
    (Agrid .* zgrid).^(p.sigma-1) ) .^ (1/denom);
prod = Agrid .* zgrid .* (kgrid.^p.alpha) .* (labor.^(1-p.alpha));
price = (Pgrid.^(p.sigma-1).*Igrid ./ prod).^(1/p.sigma);

profit = prod.*price - labor.*wgrid;

% value function iteration
tolV = 1; Vct = 0;
while tolV > 0.01
    % compute expected value of V(z',k',A',K') given (z,k) and (A,K), 
    % the transition matrices to (A', z'), as well as the (known) future K'
    % according to the regression coefficients
    expval = get_expected_value(V, K, guess, pi, p);
    
    % Now transform expectation into 5-dim format with k' on d5 and add
    helper = permute(repmat(kgrid, [1 1 1 1 p.Nk]), [1 5 3 4 2]) - ...
        (1 - p.delta) * repmat(kgrid, [1 1 1 1 p.Nk]);
    Vmat_adj = repmat(P_mugrid, [1 1 1 1 p.Nk]) .* ...
        (repmat((profit - p.phi), [1 1 1 1 p.Nk]) - helper) + ...
        p.beta * repmat(permute(expval, [1 5 3 4 2]), [1 p.Nk 1 1 1]);
    [Vnew_adj, pol_adj] = max(Vmat_adj, [], 5);
    
    % now calculate expected value when not adjusting, ie having k' = 
    % (1-delta)*k
    k_depr = (1-p.gA) * k; % is the amount of capital left if you do not want to pay fixed adjustment cost
    helper = permute(expval, [2 1 3 4]);
    expval_non = interp1(k, helper, [k(1) k_depr(2:end)], 'linear');
    expval_non = permute(expval_non, [2 1 3 4]);
    
    Vnew_non = P_mugrid .* (profit - repmat((p.delta - p.gA) * k, [p.Nz 1 p.NA p.NK])) + ...
        p.beta * expval_non;
    
    Vnew = max(max(Vnew_non, Vnew_adj), 0);
%     Vnew = max(Vnew_non, Vnew_adj - p.phi);
    
    tolV = max(abs(Vnew(:)-V(:)));
    Vct = Vct + 1;
    if mod(Vct,50)==0 
        disp([Vct, tolV])
    end
    
    V = Vnew;
end

pol_non = repmat(interp1(k, 1:p.Nk, [k(1) k_depr(2:end)],'linear'), [p.Nz 1 p.NA p.NK] );
pol_act = Vnew_adj - p.phi > Vnew_non;
pol = pol_act .* pol_adj + (1 - pol_act) .* pol_non;

pol_fl = max(floor(pol), 1);
pol_ce = min(ceil(pol), p.Nk);
pol_weight = pol - pol_fl;
polfun = (1-pol_weight) .* k(pol_fl) + pol_weight .* k(pol_ce);

adjcost = (pol_act * p.phi);
investment = (polfun - (1-p.delta) * kgrid); 
% is gross investment of an individual firm (net of adjustment costs, but including depreciation)
investment_adj = investment + adjcost; % includes adjustment cost paid
% maybe this should be called Q(k,k'(z,k,A,K))
 
toc
disp(['Vfn done, Vct= ', num2str(Vct)])

%% SIMULATION
% now simulate the economy and esp the aggregate capital stock 

% start with uniform distribution over states (z,k)
Mu = ones(p.Nz,p.Nk); Mu = Mu/sum(Mu(:));

agg_req.prod = prod;
agg_req.inv_adj = investment_adj;
agg_req.profit = profit;

agg_opt.price = price;
agg_opt.rev_prod = zgrid .* Agrid;
agg_opt.labordemand = labor;
agg_opt.num_adjusters = pol_act;


mu = ones(p.Nz, p.Nk);
mu = mu ./ sum(mu(:));
[agg_ts, final_mu] = KS_sim(mu,K,k,pi,pol,Astate,p.T,agg_req,agg_opt,p);

aggVars = zeros(p.T,4);
aggVars(:,1) = agg_ts.K;
aggVars(:,2) = agg_ts.C;
aggVars(:,3) = agg_ts.prod;
aggVars(:,4) = agg_ts.w;


% kill burn-in period:
burn_in = 500;
aggVars(1:burn_in,:) = []; % columns: K, C, I, w
agg_ts.labordemand(1:burn_in,:) = [];
agg_ts.rev_prod(1:burn_in) = [];
agg_ts.num_adjusters(1:burn_in) = [];
agg_ts.highest_k(1:burn_in,:) = [];

display(mean(aggVars));

meanlabor = zeros(1,p.NA);
for iA = 1:p.NA
    meanlabor(iA) = mean(agg_ts.labordemand(Astate(burn_in+1:end-1)==iA));
end
mean_TFPRdisp = zeros(1,p.NA);
for iA = 1:p.NA
    mean_TFPRdisp(iA) = mean(agg_ts.rev_prod(Astate(burn_in+1:end-1)==iA));
end

toc
disp('Time series of moments computed')
disp(['Mean labor demand = ', num2str(meanlabor), ' highest k: ', num2str(max(agg_ts.highest_k))])

% check if grid OK
% if max(highest_k)==length(k)
%     maxk = 2 * maxk;
%     disp(['Too much mass on top of grid -- redo Vfi with maxk = ', num2str(maxk)])
% elseif max(highest_k)<= 0.5 * length(k)
%     maxk = 2/3 * maxk;
%     disp(['Too much mass in lower half of grid -- redo Vfi with maxk = ', num2str(maxk) ])
% else % enter regression
      
%% REGRESSION
newcoeffs = zeros(size(coeffs));

Areg = Astate(burn_in+1:end-1)';
Xmat = zeros(p.T - burn_in-1, p.NA + 1);
logKvec = log(aggVars(1:end-1,1));
for iA = 1:p.NA
    Xmat(:,iA) = Areg == iA;
end
Xmat(:,end) = logKvec; % regressors: constant plus log(K)

clear Ymat
% Y vectors:
Ymat(:,1) = aggVars(2:end, 1); % is the vector of
% aggregate capital in periods _following_ periods of state i
Ymat(:,2:4) = aggVars(1:end-1, 2:4);

beta_reg = (Xmat'*Xmat)\(Xmat'*Ymat);
% beta_reg(end,3) = 0; % P = 1 fixed, independent of K
newcoeffs(:,:) = beta_reg;

coefftol = max(abs(coeffs(:)-newcoeffs(:)))
coeffct = coeffct + 1
if mod(coeffct, length(coefftol_hist)) == 0
    if mean(coefftol_hist(6:end)) > 0.95 * mean(coefftol_hist(1:5))
        p.alpha_adj = p.alpha_adj / 2;
    else
        p.alpha_adj = min(2 * p.alpha_adj, alpha_adj_init);
    end
    display(p.alpha_adj)
    coefftol_hist = coefftol_hist * 0;
else
    coefftol_hist(mod(coeffct, length(coefftol_hist))) = coefftol;
end
coeffs = p.alpha_adj*newcoeffs + (1-p.alpha_adj)*coeffs

% R_flag = 1;
% R_old = R;
% display(coeffs);
if mod(coeffct,20) == 0
    pause(5)
end

% end

end

mean_adjusters = zeros(1,p.NA);
for iA = 1:p.NA
    mean_adjusters(iA) = mean(agg_ts.num_adjusters(Astate(burn_in+1:end)==iA));
end

% IMPULSE RESPONSES!
% set up initial distribution
neutral_state = zeros(100,1) + 3;
[neutral_agg, neutral_mu] = KS_sim(final_mu,K,k,pi,pol,neutral_state,100,agg_req,agg_opt,p);

sim_dur = 20; % no. of periods to simulate following shock
num_econs = 5000; % no. of economies to simulate
response = zeros(sim_dur + 3, 4, num_econs); % K, inv_adj, 
response(1:3, 1, :) = neutral_agg.K(end);
response(1:3, 2, :) = neutral_agg.inv_adj(end);
response(1:3, 3, :) = neutral_agg.C(end);
response(1:3, 4, :) = neutral_agg.prod(end);

Acdfs = cumsum(pi.A,2);
sim_state = zeros(sim_dur, num_econs);
sim_state(1,:) = 1;
counter = 0;
for ni = 1:num_econs
    seed = ni; rng(seed);
    sim_shocks = rand(sim_dur-1,1);
    for t = 1:sim_dur-1
        sim_state(t+1,ni) = find(sim_shocks(t) <= Acdfs(sim_state(t, ni),:), 1, 'first') ;
    end
    ind = find(sim_state(2:end, ni) >=3, 1, 'first');
    sim_state(1+ind:end,ni) = 3;
    
    [eco_agg, ~] = KS_sim(neutral_mu,K,k,pi,pol,sim_state,sim_dur,agg_req,agg_opt,p);
    response(4:end, 1, ni) = eco_agg.K;
    response(4:end, 2, ni) = eco_agg.inv_adj;
    response(4:end, 3, ni) = eco_agg.C;
    response(4:end, 4, ni) = eco_agg.prod;
    
    % make cool looking counter
    if ni >= (counter + 1) * num_econs / 20
        fprintf('.')
        counter = counter + 1;
    end
    
end
fprintf('\n')





