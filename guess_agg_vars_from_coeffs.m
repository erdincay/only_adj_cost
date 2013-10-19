function guess = guess_agg_vars_from_coeffs(p, coeffs, K, pi)
maxK = max(K);
minK = min(K);

guess.Kprime = zeros(p.NA,p.NK);
guess.C = guess.Kprime;
guess.I = guess.Kprime;
guess.w = guess.Kprime;
for i = 1:p.NA
    guess.Kprime(i,:) = coeffs(1,1,i) + log(K) * coeffs(2,1,i);
    guess.C(i,:)      = coeffs(1,2,i) + log(K) * coeffs(2,2,i);
    guess.P(i,:)      = coeffs(1,3,i) + log(K) * coeffs(2,3,i);
    guess.I(i,:)      = coeffs(1,4,i) + log(K) * coeffs(2,4,i);
    guess.w(i,:)      = coeffs(1,5,i) + log(K) * coeffs(2,5,i);
end

lownumber = 0.000001;
if any(guess.Kprime(:) < 0)
    guess.Kprime(guess.Kprime < minK) = minK;
    disp(['Negative values for K'' adjusted to ', num2str(minK)])
end
if any(guess.Kprime(:) > maxK)
    guess.Kprime(guess.Kprime > maxK) = maxK;
    disp(['Too high values for K'' adjusted to ', num2str(maxK)])
end
if any(guess.P(:) < 0)
    guess.P(guess.P < 0) = lownumber;
    disp(['Negative values for P adjusted to ', num2str(lownumber)])
end
if any(guess.C(:) < 0)
    guess.C(guess.C < 0) = lownumber;
    disp(['Negative values for C adjusted to ', num2str(lownumber)])
end
if any(guess.I(:) < 0)
    guess.I(guess.I < 0) = lownumber;
    disp(['Negative values for I adjusted to ', num2str(lownumber)])
end
if any(guess.w(:) < 0)
    guess.w(guess.w < 0) = lownumber;
    disp(['Negative values for w adjusted to ', num2str(lownumber)])
end

guess.R = zeros(p.NA,p.NK);
mUC = guess.C .^ (-p.psi); % marginal utility of consumption
expmU = zeros(p.NA, p.NK); % expected marginal utility
for iA = 1:p.NA
    for iK = 1:p.NK
        futureC = interp1(K, permute(guess.C, [2 1]), guess.Kprime(iA,iK),'linear', 'extrap')';
        expmU(iA,iK) = pi.A(iA,:) * (futureC .^ (-p.psi));
        guess.R(iA,iK) = mUC(iA,iK) / (p.beta * expmU(iA,iK));
    end
end