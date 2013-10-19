function expval = get_expected_value(V, K, guess, pi, p)

expval = zeros(p.Nz,p.Nk,p.NA,p.NK);
    
for toK = 1:p.NK
    for fromz = 1:p.Nz
        expvalz = zeros(p.Nk,p.NA); % is expected value of being in 'fromz' today

        for toz = 1:p.Nz
            expvalA = pi.A*permute(V(toz,:,:,toK), [3 2 1 4]);
            expvalz = expvalz + pi.z(fromz,toz)*expvalA';
        end

        expval(fromz,:,:,toK) = expvalz;
    end
end
% this is the expected value of being in state (z,A) conditional 
% on going to (k',K') --> ie (z,k',A,K'). Use coefficients from regressing K' 
% on K and interpolate. Loop over all A,K to find K' and then
% interpolate in (z,k,A,:).
helper = permute(expval, [4 1 2 3]); % move K into 1st dim for interp1 --> (K',z,k,A)
result = zeros(size(helper));
for fromK = 1:p.NK
    for fromA = 1:p.NA
        result(fromK,:,:,fromA) = interp1(K, helper(:,:,:,fromA), guess.Kprime(fromA,fromK), 'linear','extrap' );
    end
end
expval = permute(result, [2 3 4 1]); % move K back to 4th dim - this is
% the expectation of being in (z,A,K) conditional on going to k' -->
% (z,k',A,K)
