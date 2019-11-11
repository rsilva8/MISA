function combinatorial_optim(O, varargin)
% Solve combinatorial optimization

if ~isempty(varargin)
    w0 = O.greedysearch(varargin{1});
else
    w0 = O.greedysearch();
end


% w0 = O.sub_perm_analysis(w0);

O.objective(w0);

end