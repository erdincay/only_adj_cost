function [coeffs, mean_adjusters, mean_TFPRdisp] = adj_cost_fun(p)

% PARAMETERS
% p.beta = 0.98; % discount factor
% p.delta = 0.16; % depreciation
% p.alpha = 0.33; % production function: y = A*z * k^alpha * l^(1-alpha)
% p.rhoA = 0.9; % persistence in A
% p.sigA = 0.03; % stddev in shocks to A
% p.gA = 0.01; % growth rate in A (use to ensure eventual adjustment)
% p.rhoz = 0.7; % persistence in z
% p.sigz = 0.6; % stddev in shocks to z
% p.sigma = 2; % elasticity of substitution between goods
% p.phi = 0.5; % fixed adjustment cost (0 = frictionless)
% p.psi = 1.05; % CRRA coefficient on household utility
% 
% p.sigexp = (p.sigma-1)/p.sigma;
% 
% % program
% p.Nk = 100; p.maxk = 1.3; p.mink = 0; % grid capital
maxk = p.maxk; mink = p.mink;
% p.NA = 5;
% p.Nz = 7;
% p.NK = 50; 
% p.K_range_lo = 0.8; p.K_range_hi = 1.3;
maxK = p.K_range_hi*(maxk-mink); minK = p.K_range_lo*(maxk-mink);
% p.alpha_adj = 0.4; % weight on new coefficients
% p.T = 1700; % length of simulation

% flags
flags.new_calc = 1;
flags.new_shocks = 0;
flags.matrices_to_update_Mu = 0;

% Set up grids
k = mink:(maxk-mink)/(p.Nk):maxk; k(1)=[];
K = minK:(maxK-minK)/(p.NK):maxK; K(1)=[];
[A, pi.A] = TAUCHEN(p.NA, p.rho_A, p.sig_A, p.A_tauchen); A = exp(A');
[z, pi.z] = TAUCHEN(p.Nz, p.rho_z, p.sig_z, p.z_tauchen); z = exp(z');

% construct vector of aggregate shocks:
Astate = aggregate_shocks(pi, p, flags);

% get guesses from economy without aggregate uncertainty
noA = adj_cost_noA(p);
display('Initial values fetched');

% Set up grids not relying on k
% V = zeros(p.Nz, p.Nk, p.NA, p.NK);
V = repmat(noA.V, [1 1 p.NA p.NK]);
maxk = p.K_range_hi * A(end) * z(end) * noA.k(find(sum(noA.Mu),1,'last'));
% mink = p.K_range_lo * z(1) * noA.k(find(sum(noA.Mu),1,'first'));

zgrid = repmat(z', [1 p.Nk p.NA p.NK]);
kgrid = repmat(k, [p.Nz 1 p.NA p.NK]);
Agrid = repmat(permute(A, [1 3 2]), [p.Nz p.Nk 1 p.NK]);

% Initial values for coefficients
coeffs = zeros(2, 5, p.NA); N_coeff = 5; % K, C, P, I, w
coeff_hist = zeros(2, N_coeff, p.NA, 10);
coeffs(1,1,:) = noA.K_agg; % agg_K;
coeffs(1,2,:) = noA.C_agg; % C;
% coeffs(1,3,:) = noA.P; % P;
coeffs(1,3,:) = 1; % P;
coeffs(1,4,:) = noA.I_agg; % I;
coeffs(1,5,:) = noA.wguess; % w;

R_old = zeros(p.NA,p.NK) + 1 / p.beta;
R_flag = 1;
alpha_R = 0.05;

% Start coefficients loop
coefftol = 1; coeffct = 0;
while coefftol > 0.005 && coeffct < 50

% Set up grids relying on k
k = mink:(maxk-mink)/(p.Nk):maxk; k(1)=[];
maxK = p.K_range_hi * noA.K_agg; minK = p.K_range_lo*noA.K_agg;
K = minK:(maxK-minK)/(p.NK):maxK; K(1)=[];

kgrid = repmat(k, [p.Nz 1 p.NA p.NK]);

% get guesses from coefficients
guess = guess_agg_vars_from_coeffs(p, coeffs, K, pi);

R = zeros(p.NA, p.NK);
mUC = guess.C .^ (-p.psi); % marginal utility of consumption
expmU = zeros(p.NA, p.NK); % set up expected marginal utility
for iA = 1:p.NA
    for iK = 1:p.NK
        futureC = interp1(K, permute(guess.C, [2 1]), guess.Kprime(iA,iK),'linear', 'extrap')';
        expmU(iA,iK) = pi.A(iA,:) * (futureC .^ (-p.psi));
        R(iA,iK) = mUC(iA,iK) / (p.beta * expmU(iA,iK));
    end
end
if R_flag % && coefftol < 0.1
    R_tol = max(abs(R(:) - R_old(:)))
    R = alpha_R * R + (1 - alpha_R) * R_old
    if any(R(:) < 1)
        disp('adjusting R<1')
%         R = max(R,1);
    end
else
    R = R_old;
end
R_flag = 0;

Pgrid = repmat(permute(guess.P, [3 4 1 2]), [p.Nz p.Nk 1 1]);
Igrid = repmat(permute(guess.I, [3 4 1 2]), [p.Nz p.Nk 1 1]);
wgrid = repmat(permute(guess.w, [3 4 1 2]), [p.Nz p.Nk 1 1]);
Rgrid = repmat(permute(R, [3 4 1 2]), [p.Nz p.Nk 1 1]);

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
    Vmat_adj = repmat(profit, [1 1 1 1 p.Nk]) - helper + ...
        repmat(permute(expval ./ Rgrid, [1 5 3 4 2]), [1 p.Nk 1 1 1]);
    [Vnew_adj, pol_adj] = max(Vmat_adj, [], 5);
    
    % now calculate expected value when not adjusting, ie having k' = 
    % (1-delta)*k
    k_depr = (1 - p.gA) * k; % is the amount of capital left if you do not want to pay fixed adjustment cost
    helper = permute(expval, [2 1 3 4]);
    expval_non = interp1(k, helper, k_depr, 'linear', 'extrap');
    expval_non = permute(expval_non, [2 1 3 4]);
    
    Vnew_non = profit - repmat((p.delta - p.gA) * k, [p.Nz 1 p.NA p.NK]) + ...
        expval_non ./ Rgrid;
    
    Vnew = max(max(Vnew_non, Vnew_adj - p.phi), 0);
%     Vnew = max(Vnew_non, Vnew_adj - p.phi);
    
    tolV = max(abs(Vnew(:)-V(:)));
    Vct = Vct + 1;
    if mod(Vct,50)==0 
        disp([Vct, tolV])
    end
    
    V = Vnew;
    if Vct > 1000
        error('Vfn did not converge');
    end
end

pol_non = repmat(interp1(k, 1:p.Nk, k_depr,'linear','extrap'), [p.Nz 1 p.NA p.NK] );
pol_act = Vnew_adj - p.phi > Vnew_non;
pol = pol_act.*pol_adj + (1-pol_act).*pol_non;

pol_fl = max(floor(pol), 1);
pol_ce = min(ceil(pol), p.Nk);
pol_weight = pol - pol_fl;
polfun = (1-pol_weight) .* k(pol_fl) + pol_weight .* k(pol_ce);

adjcost = (pol_act.*p.phi);
investment = (polfun - (1-p.delta) * kgrid); 
% is gross investment of an individual firm (net of adjustment costs, but including depreciation)
investment_adj = investment + adjcost; % includes adjustment cost paid
% maybe this should be called Q(k,k'(z,k,A,K))
 
disp(['Vfn done, Vct= ', num2str(Vct)])

%% SIMULATION
% now simulate the economy and esp the aggregate capital stock 

% start with uniform distribution over states (z,k)
Mu = ones(p.Nz,p.Nk); Mu = Mu/sum(Mu(:));

% make list to keep track of aggregate variables
aggVars = zeros(p.T,4); % columns: K, C, I, w
labordemand = zeros(p.T,1);
num_adjusters = zeros(p.T,1);
highest_k = zeros(p.T,1);% keep track of the highest observed amount of 
% capital held. If it is maxk too frequently, agents are likely to be 
% constrained by the choice of the grid --> need to change the grid

rev_prod = zeros(p.T,1);
Mu_hist = zeros(p.Nz, p.Nk, p.T);

for t = 1:p.T
    % compute aggregate variables
    K_t = sum(Mu) * k'; % aggregate capital
    Kt_ind = find(K > K_t, 1, 'first'); % index of agg. capital in K
    if isempty(Kt_ind)
        Kt_ind = length(K)+1;
    end
    aggVars(t,1) = K_t;
    
    Kt_ind_hi = min(Kt_ind, length(K));
    Kt_ind_lo = max(1, Kt_ind-1);
    if Kt_ind_hi > Kt_ind_lo
        weight_hi = (K_t - K(Kt_ind_lo)) / (K(Kt_ind_hi) - K(Kt_ind_lo)) ; 
    else
        weight_hi = 1;
    end
    % is the weight on the upper of the two grid points in K that Kt falls 
    % in between of. 
    
    % Production Y
    prod_dist = (1 - weight_hi) * prod(:,:,Astate(t),Kt_ind_lo) + ...
        weight_hi * prod(:,:,Astate(t),Kt_ind_hi);
    Y_t = (sum( (prod_dist(:).^ p.sigexp) .* Mu(:) )) ^ (1/p.sigexp);
    
    % Investment and adjustment cost
    inv_adj_dist = (1 - weight_hi) * investment_adj(:,:,Astate(t),Kt_ind_lo) + ...
        weight_hi * investment_adj(:,:,Astate(t),Kt_ind_hi);
    Inv_t = sum( inv_adj_dist(:) .* Mu(:) ); % is gross investment plus
    % adjustment cost: k' - (1-delta)k (+ phi if nec) 
    
    % Consumption C
    C_t = Y_t - Inv_t;  
    
%     % Price index P - don't need if P=1 fixed
    price_dist = (1 - weight_hi) * price(:,:,Astate(t),Kt_ind_lo) + ...
        weight_hi * price(:,:,Astate(t),Kt_ind_hi);
%     P_t = sum( (price_dist(:).^(1-p.sigma)).* Mu(:) ).^(1/(1-p.sigma));
    P_t = 1;
    
    % expenditures for intermediates I
    I_t = Y_t * P_t;
    
    % wage w
    % agg wage bill (divided by 1 unit of labor) = agg income - agg profits
    profit_dist = (1 - weight_hi) * profit(:,:,Astate(t),Kt_ind_lo) + ...
        weight_hi * profit(:,:,Astate(t),Kt_ind_hi);
    
    profit_t = sum( profit_dist(:) .* Mu(:) ) - Inv_t;
    
    w_t = C_t - profit_t;

    aggVars(t,1) = K_t;
    aggVars(t,2) = C_t;
    aggVars(t,3) = P_t;
    aggVars(t,4) = I_t;
    aggVars(t,5) = w_t;
    
    % revenue productivity: coefficient of variation
    rev_prod_dist = price_dist .* zgrid(:,:,1,1) .* A(Astate(t));
    [rev_histo, order] = sort(rev_prod_dist(:));    
    rev_histo_weights = Mu(order);
    rev_prod(t) = coeff_var(rev_histo, rev_histo_weights);
    
    % number of adjusters:
    adj_dist = (1 - weight_hi) * pol_act(:,:,Astate(t),Kt_ind_lo) + ...
        weight_hi * pol_act(:,:,Astate(t),Kt_ind_hi);
    num_adjusters(t) = sum(adj_dist(:) .* Mu(:));
  
%     display(aggVars(t,:));
    
    % Labor demand
    labordist = (1-weight_hi)*labor(:,:,Astate(t),Kt_ind_lo) +...
        weight_hi * labor(:,:,Astate(t),Kt_ind_hi);
    labordemand(t) = sum( labordist(:) .* Mu(:) )  ;
    
    % save Mu in history
    Mu_hist(:,:,t) = Mu;
    highest_k(t) = find(sum(Mu),1,'last');
    
    % update distribution
    MuNew = update_dist(Mu, pi, pol, Astate(t), Kt_ind_hi, Kt_ind_lo, weight_hi);
    
    Mu = MuNew;
end

% kill burn-in period:
burn_in = 700;
aggVars(1:burn_in,:) = []; % columns: K, C, I, w
labordemand(1:burn_in,:) = [];
rev_prod(1:burn_in) = [];
num_adjusters(1:burn_in) = [];
highest_k(1:burn_in,:) = [];

display(mean(aggVars));

meanlabor = zeros(1,p.NA);
for i = 1:p.NA
    meanlabor(i) = mean(labordemand(Astate(burn_in+1:end-1)==i));
end
mean_TFPRdisp = zeros(1,p.NA);
for i = 1:p.NA
    mean_TFPRdisp(i) = mean(rev_prod(Astate(burn_in+1:end-1)==i));
end
mean_adjusters = zeros(1,p.NA);
for iA = 1:p.NA
    mean_adjusters(iA) = mean(num_adjusters(Astate(burn_in+1:end)==iA));
end

disp('Time series of moments computed')
disp(['Mean labor demand = ', num2str(meanlabor), ' highest k: ', num2str(max(highest_k))])
disp(['Mean TFPR disp = ', num2str(mean_TFPRdisp)])
disp(['Mean adjusters = ', num2str(mean_adjusters)])

% check if grid OK
if max(highest_k)==length(k)
    maxk = 2 * maxk;
    disp(['Too much mass on top of grid -- redo Vfi with maxk = ', num2str(maxk)])
elseif max(highest_k)<= 0.5 * length(k)
    maxk = 2/3 * maxk;
    disp(['Too much mass in lower half of grid -- redo Vfi with maxk = ', num2str(maxk) ])
else % enter regression
      
%% REGRESSION
newcoeffs = zeros(size(coeffs));
for i = 1:p.NA
    Xvec = aggVars(Astate(burn_in+1:end-1)==i, 1);
    Xmat = [ones(size(Xvec)) log(Xvec)]; % regressors: constant plus log(K)
    
    clear Ymat
    % Y vectors:
    Ymat(:,1) = aggVars(logical( [0 Astate(burn_in+1:end-1)==i] ) , 1)'; % is the vector of 
    % aggregate capital in periods _following_ periods of state i
    Ymat(:,2) = aggVars(Astate(burn_in+1:end-1)==i, 2)';
    Ymat(:,3) = aggVars(Astate(burn_in+1:end-1)==i, 3)';
    Ymat(:,4) = aggVars(Astate(burn_in+1:end-1)==i, 4)';
    Ymat(:,5) = aggVars(Astate(burn_in+1:end-1)==i, 5)';
    
    beta_reg = (Xmat'*Xmat)\(Xmat'*Ymat);
    newcoeffs(:,:,i) = beta_reg;
    
    % if P=1 fixed:
    newcoeffs(:,3,i) = [1 0]';
    
%     disp(corr(Xmat * beta_reg,Ymat))
end

coefftol = max(abs(coeffs(:)-newcoeffs(:)))
coeffct = coeffct + 1
coeffs = p.alpha_adj*newcoeffs + (1-p.alpha_adj)*coeffs;

R_flag = 1;
R_old = R;
% display(coeffs);
if mod(coeffct,20) == 0
    pause(5)
end

end

end


disp('Done')
