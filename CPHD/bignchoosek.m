function [out] = bignchoosek(n,k)
% Ramanujan's approximation for n choose k
% https://stackoverflow.com/questions/40527010/r-how-can-i-calculate-large-numbers-in-n-choose-k
%out = exp(ramanujan(n) - ramanujan(k) - ramanujan(n-k));

out = exp( gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1) );
end

function v = ramanujan(n)
v = n*log(n) - n + log(n*(1 + 4*n*(1+2*n)))/6 + log(pi)/2;
end