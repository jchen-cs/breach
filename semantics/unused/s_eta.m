% The eta function, which is an optimistic integrator typically used in
% UNTIL for rho-
% That is, it takes the "better" score of its two inputs.

% For vectors of robustness values v1 and v2 for STL formulae phi1 and
% phi2, compute eta(v1, v2) according to the name of the semantics 
% requested.

% eta has the same soundness conditions as beta.

% When used in rho-, eta is usually negated (ie -eta(-v1, -v2) ) so it's
% treated as a pessimistic integrator.

function res = s_eta(v1, v2, semantics)
    switch semantics
        case 'max-breach'
            res = max(v1, v2);
        case 'const-breach'
            res = max(v1, v2);
        case 'plus-breach'
            res = max(v1, v2);
        case 'telex'
            res = max(v1, v2);
        case 'belta'
            res = max(v1, v2);
        case 'agm-product'
            res = max(v1, v2);
        case 'sum-product'
            res = v1 + v2;
        case 'sum-min'
            res = v1 + v2;
        case 'max-product'
            res = max(v1, v2);
        case 'minonly'
            res = min(v1, v2);
        case 'smoothrect'
            res = max(v1, v2);
        otherwise
            error('Unknown semantics for eta!');
    end
end