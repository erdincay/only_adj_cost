function [y,yprob] =TAUCHEN(ny,lambda,sigma,m)
% TAUCHEN Tauchen's algorithm (1986)
%       A Markov chain whose sample paths approximate those of the AR(1) process 

%       y(t+1) =lambda * y(t) + eps(t+1)

%       where eps are normal with std. deviation sigma

%     ny                scalar, number of points in y-grid
%     lambda            scalar
%     sigma             scalar, std. dev. of epsilons
%     m                 max +- std. devs.     
%     y                 1 * ny vector, grid for y
%     yprob             ny * ny matrix, transition probabilities

%The values of the following parameters are chosen in the light of Tauchen's paper

% sigma=.011;
% m=0.011;
% ny=4;
% lambda=0.95;

y     = zeros(ny,1);
yprob = zeros(ny,ny);

y(ny) = m * sqrt(sigma^2 / (1 - lambda^2));
y(1)  = -y(ny);
ystep = (y(ny) - y(1)) / (ny - 1);

for i=2:(ny-1)
    y(i) = y(1) + ystep * (i - 1);
end 

for j = 1:ny
    for k = 1:ny
        if k == 1
            yprob(j,k) = normcdf((y(1)- lambda * y(j) + ystep / 2) / sigma);
        elseif k == ny
            yprob(j,k) = 1 - normcdf((y(ny)- lambda * y(j) - ystep / 2) / sigma);
        else
            yprob(j,k) = normcdf((y(k)- lambda * y(j) + ystep / 2) / sigma) - ...
                         normcdf((y(k)- lambda * y(j) - ystep / 2) / sigma);
        end
    end
end

%   Remark: function c = normcdf(x);       where   c = 0.5 * erfc(-x/sqrt(2));   See help menu for "normcdf"
   
y;
yprob;