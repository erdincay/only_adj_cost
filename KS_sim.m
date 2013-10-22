function [ new_aggs, mu ] = KS_sim( mu, K, k, pi, pol, Astate, T, agg_req, agg_opt, p)
%KS_sim Krussel-Smith simulation of the economy
%   

new_aggs.K = zeros(T,1);
new_aggs.C = zeros(T,1);
new_aggs.w = zeros(T,1);
new_aggs.highest_k = zeros(T,1);

agg_names = [fieldnames(agg_req)', fieldnames(agg_opt)'];
for agg_name = agg_names
    new_aggs.(char(agg_name)) = zeros(T, 1);
end

for t = 1:T
    K_t = sum(mu) * k'; % aggregate capital
    Kt_ind = find(K > K_t, 1, 'first'); % index of agg. capital in K
    if isempty(Kt_ind)
        Kt_ind = length(K)+1;
    end
    new_aggs.K(t) = K_t;
    
    Kt_ind_hi = min(Kt_ind, length(K));
    Kt_ind_lo = max(1, Kt_ind-1);
    if Kt_ind_hi > Kt_ind_lo
        weight_hi = (K_t - K(Kt_ind_lo)) / (K(Kt_ind_hi) - K(Kt_ind_lo)) ; 
    else
        weight_hi = 1;
    end
    
    agg_state = [Astate(t), Kt_ind_lo, Kt_ind_hi, weight_hi];
    
    % handle required and extra variables
    prod_dist = (1 - weight_hi) * agg_req.prod(:,:,Astate(t),Kt_ind_lo) + ...
        weight_hi * agg_req.prod(:,:,Astate(t),Kt_ind_hi);
    new_aggs.prod(t) = (sum( (prod_dist(:).^ p.sigexp) .* mu(:) )) ^ (1/p.sigexp);
    new_aggs.inv_adj(t) = compute_aggregate(mu, agg_req.inv_adj, agg_state); % gross investment including adj costs!
    new_aggs.profit(t) = compute_aggregate(mu, agg_req.profit, agg_state) - new_aggs.inv_adj(t);
    new_aggs.C(t) = new_aggs.prod(t) - new_aggs.inv_adj(t);
    new_aggs.w(t) = new_aggs.C(t) - new_aggs.profit(t);
    
    % optional variables
    for agg_name = fieldnames(agg_opt)'
        switch char(agg_name)
            case 'rev_prod' % requires agg_opt.price
                price_dist = (1 - weight_hi) * agg_opt.price(:,:,Astate(t),Kt_ind_lo) + ...
                    weight_hi * agg_opt.price(:,:,Astate(t),Kt_ind_hi);
                rev_prod_dist = price_dist .* agg_opt.rev_prod(:,:,Astate(t), Kt_ind_lo); % agg_opt.rev_prod only Agrid * zgrid...
                [rev_histo, order] = sort(rev_prod_dist(:));
                rev_histo_weights = mu(order);
                new_aggs.rev_prod(t) = coeff_var(rev_histo, rev_histo_weights);
            otherwise
                new_aggs.(char(agg_name))(t) = ...
                    compute_aggregate(mu, agg_opt.(char(agg_name)), agg_state);
        end
    end
    
    new_aggs.highest_k(t) = find(sum(mu), 1, 'last');
    
    mu_new = update_dist(mu, pi, pol, agg_state);
    mu = mu_new;
    
end



end

