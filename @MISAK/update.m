function update(O,S,M,b,l,e)

O.S = S;
O.M = M;

O.d = full(sum([O.S{O.M}],2));  % Dimensionality of each subspace
O.nes = O.d~=0;             % Non-empty subspace indexes
O.d_k = cellfun(@(s) full(sum(s,2)), O.S,'Un',0);

O.K = sum(O.nes);

O.beta = b;
O.lambda = l;
O.eta = e;
O.nu = (2*O.eta + O.d - 2)./(2*O.beta);
O.a = (O.lambda.^(-1./(O.beta)) .* gamma(O.nu + 1./O.beta)) ./ ...
    (O.d .* gamma(O.nu));

end