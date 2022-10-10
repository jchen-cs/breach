% The beta function, which an optimistic integrator typically used for 
% AND in rho-
% That is, it takes the "better" score of its two inputs.

% For vectors of robustness values v1 and v2 for STL formulae phi1 and
% phi2, compute beta(v1, v2) according to the name of the semantics 
% requested.

% beta is sound when, for x1 and x2 >= 0, beta(x1, x2) > 0 only if 
% x1 > 0 or x2 > 0. 0 otherwise.

% When used in rho-, beta is usually negated (ie -beta(-v1, -v2) ) so it's
% treated as a pessimistic integrator.

function res = s_beta(v1, v2, semantics)
    switch semantics
        case 'max-breach'
            res = max(v1, v2);
        case 'const-breach'
            res = max(v1, v2);
        case 'plus-breach'
            res = koenAndPlus(v1, v2);
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
        case 'smooth1'
            res = v1 + v2;
        case 'smooth2'
            %res = log(exp(v1) + exp(v2)) / 2;
            res = v1 + v2;
        otherwise
            error('Unknown semantics for beta!');
    end
end