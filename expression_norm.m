function nEXP = expression_norm(T,logEXP,I,E)

if (nargin < 4)
    E = 2;
end

n = size(logEXP,1);
ut = unique(T);

M = zeros(size(T));
for i = 1:max(size(ut))
    x = logEXP(I,T==ut(i));
    M(T==ut(i)) = mean(x(:));%log2(mean(2.^x(:)));
end

%M = M - mean(M);
M = M - E;
nEXP = logEXP - repmat(M,n,1);
