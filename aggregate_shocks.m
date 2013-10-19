function Astate = aggregate_shocks(pi, p, flags)

if flags.new_shocks == 1 % draw new random numbers
    % construct vector of aggregate shocks:
    randostring = rand(1,p.T);
    save randostring randostring
    disp('Drawn standardnormal random numbers')
else  % load Astate vector
    load randostring
end
Astate = zeros(1,p.T); 
Astate(1) = 1; % start in the lowest state
Acdfs = cumsum(pi.A,2);
for i = 1:p.T-1
    Astate(i+1) = find(randostring(i) <= Acdfs(Astate(i),:), 1, 'first') ;
end