function CV = coeff_var(values, weights)
% weights have to sum to 1! 
values = values(weights>0);
weights = weights(weights>0);
CV = sqrt(var(values, weights)) / (values' * weights);