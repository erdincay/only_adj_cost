function agg = compute_aggregate( mu, aggregate, agg_state )
%compute_aggregates computes values of aggregate variables
%   inputs: 
% mu - a distribution over idiosyncratic states
% aggregates - a map from a firm state into a value which is then 
% aggregated according to aggregate state and mu
% agg_state - information about aggregate technology A and capital K
%   outputs:
% agg_values - the levels of the aggregate values

Astate = agg_state(1);
Kt_ind_lo = agg_state(2);
Kt_ind_hi = agg_state(3);
weight_hi = agg_state(4);

agg_dist = (1 - weight_hi) * aggregate(:,:,Astate,Kt_ind_lo) + ...
    weight_hi * aggregate(:,:,Astate,Kt_ind_hi);
agg = agg_dist(:)' * mu(:);

end

