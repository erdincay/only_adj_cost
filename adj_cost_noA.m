function noA = adj_cost_noA(p)

% PARAMETERS
% p.beta = 0.98; % discount factor
% p.delta = 0.15; % depreciation
% p.alpha = 0.33; % production function: y = k^alpha*(zAl)^(1-alpha)
% p.rhoz = 0.7; % persistence in z
% p.sigz = 0.4; % stddev in shocks to z
% p.sigma = 2; % elasticity of substitution between goods
% p.phi = 0.01; % fixed adjustment cost (0 = frictionless)
% p.psi = 2; % CRRA coefficient on household utility

p.sigexp = (p.sigma-1)/p.sigma;

% program
% p.Nk = 50; p.maxk = 1; p.mink = 0; 
maxk = p.maxk; mink = p.mink; % grid capital
% p.Nz = 4;
P = 1;

% Set up grids not relying on maxk
V = zeros(p.Nz, p.Nk);
[z, pi.z] = TAUCHEN(p.Nz, p.rho_z, p.sig_z, p.z_tauchen); z = exp(z');
[pi.z, ~, ~] = TAUCHEN_Bloom(p.Nz, z(1), z(end), 1, p.sig_z, p.rho_z, 2, 0,1);
Pguess = 1;
Iguess = 1;
wguess = 0.4;

while true % loop for maxk

% set up grids relying on maxk
k = mink:(maxk-mink)/(p.Nk):maxk; k(1)=[];
zgrid = repmat(z', [1 p.Nk]);
kgrid = repmat(k, [p.Nz 1]);

prices_tol=1; prices_ct = 1;
while prices_tol > 0.001 && prices_ct <= 100

% calculate in-period returns for all states
denom = p.alpha*(p.sigma-1)+1;
labor = ( ((p.sigma-1)/p.sigma*(1-p.alpha)./wguess).^(p.sigma) .* (Iguess.*P.^(p.sigma-1)) .* ...
    kgrid.^(p.alpha*(p.sigma-1)) .* ...
    (zgrid).^(p.sigma-1) ) .^ (1/denom);
prod = zgrid .* (kgrid.^p.alpha) .* (labor.^(1-p.alpha));
price = (P.^(p.sigma-1).*Iguess ./ prod).^(1/p.sigma);

profit = prod.*price - labor.*wguess;

% value function iteration
tolV = 1; Vct = 0;
while tolV > 0.01
    % expected value given z and k'
    expval = pi.z*V;
    
    % first case: adjustment
    ret_adj = repmat(profit, [1 1 p.Nk]) - p.phi - ...
        repmat(permute(k' * ones(1,p.Nk) - (1 - p.delta) * ones(p.Nk,1)*k, [3 2 1]), [p.Nz 1 1]); 
    Vmat_adj = ret_adj + p.beta*repmat(permute(expval, [1 3 2]), [1 p.Nk 1]);
    [Vnew_adj, pol_adj] = max(Vmat_adj, [], 3);
    
    % now calculate expected value when not adjusting, ie having k' = 
    % (1-delta)*k
    k_depr = (1-p.gA) * k; % is the amount of capital left if you do not want to pay fixed adjustment cost
    helper = permute(expval, [2 1]);
    Vnew_non = profit - repmat( (p.delta - p.gA) * k, [p.Nz 1] ) + ...
        permute(p.beta*interp1(k, helper, k_depr, 'linear', 'extrap'), [2 1]);
    
    Vnew = max(Vnew_non, Vnew_adj);
    
    tolV = max(abs(Vnew(:)-V(:)));
    Vct = Vct+1;
%     if mod(Vct,50)==0 
%         disp([Vct, tolV])
%     end
    
    V = Vnew;
end

pol_non = repmat(interp1(k, 1:p.Nk, k_depr,'linear','extrap'), [p.Nz 1] );
pol_act = Vnew_adj > Vnew_non;
pol = pol_act.*pol_adj + (1-pol_act).*pol_non;

pol_fl = max(floor(pol), 1);
pol_ce = ceil(pol);
pol_weight = pol - pol_fl; % is weight on pol_ce
polfun = (1-pol_weight) .* k(pol_fl) + pol_weight .* k(pol_ce);

adjcost = (pol_act.*p.phi);
investment = (polfun - (1-p.delta) * kgrid); 
% is gross investment of an individual firm (net of adjustment costs, but including depreciation)
investment_adj = investment + adjcost; % includes adjustment cost paid
% maybe this should be called Q(k,k'(z,k,A,K))
 

%% TRANSITIONS
Tmats = zeros(p.Nk,p.Nk,p.Nz,p.Nz);

for fromk = 1:p.Nk
    for fromz = 1:p.Nz
        Tmats(fromk, pol_fl(fromz,fromk), fromz, :) = (1-pol_weight(fromz,fromk)) * pi.z(fromz,:);
        Tmats(fromk, pol_ce(fromz,fromk), fromz, :) = Tmats(fromk,pol_ce(fromz,fromk),fromz,:) + ...
            permute(pol_weight(fromz,fromk) * pi.z(fromz,:), [1 4 3 2]) ;
    end
end

Mu = ones(p.Nz,p.Nk); Mu = Mu/sum(Mu(:));
% k is along 2nd dimension here.
mutol = 1;
while mutol > 0.00001
    MuNew = zeros(size(Mu));
    for fromz = 1:p.Nz
        for toz = 1:p.Nz
            MuNew(toz,:) = MuNew(toz,:) + Mu(fromz,:)*Tmats(:,:,fromz,toz);
        end
    end
    mutol = max(abs(Mu(:)-MuNew(:)));
    Mu = MuNew;
end

% statistics
labor_agg = sum( labor(:) .* Mu(:) );  % labor demand
K_agg = sum( kgrid(:) .* Mu(:) ); 
inv_adj_agg = sum( investment_adj(:) .* Mu(:) );
prod_agg = (sum( (prod(:).^ p.sigexp) .* Mu(:) )) ^ (1/p.sigexp);
P_index = sum( (price(:).^(1-p.sigma)).* Mu(:) ).^(1/(1-p.sigma));
I_agg = prod_agg * P_index;
profit_net_agg = sum( profit(:) .* Mu(:) ) - inv_adj_agg; % intra-period minus inv_adj
C_agg = prod_agg - inv_adj_agg;
% disp([labor_agg K_agg inv_adj_agg prod_agg C_agg]);

% adjust guesses
w_new = (1 + 0.2*(labor_agg-1)) * wguess;
P_new = 0.8*Pguess + 0.2*P_index;
I_new = 0.8*Iguess + 0.2*I_agg;

prices_ct  = prices_ct + 1;
prices_tol = max(abs([(w_new-wguess)/wguess; (I_new-Iguess)/Iguess; (P_index-Pguess)/Pguess; labor_agg-1]));

Pguess = P_new;
Iguess = I_new;
wguess = w_new;
% disp([Pguess Iguess wguess]);

end

maxk_ind = find(sum(Mu, 1), 1, 'last');
if maxk_ind == p.Nk 
    maxk_new = 2*maxk;
elseif maxk_ind <= 0.6*p.Nk
    maxk_new = 0.75* maxk;
else
    break
end

maxk = maxk_new;

end

noA.V = V;
noA.pol = pol;
noA.polfun = polfun;
noA.Mu = Mu;
noA.K_agg = K_agg;
noA.C_agg = C_agg;
noA.P = P_index;
noA.Y = prod_agg;
noA.I_agg = I_new;
noA.inv_adj_agg = inv_adj_agg;
noA.profit_net_agg = profit_net_agg;
noA.wguess = wguess;
noA.labor = labor_agg;
noA.maxk = maxk;
noA.prices_tol = prices_tol;
noA.k = k;

display(prices_tol)
display(noA)