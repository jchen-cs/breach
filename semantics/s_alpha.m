% The alpha function, which is the pessimistic integrator typically used
% for AND in rho+
% That is, it takes the "worse" score of its two inputs.

% For vectors of robustness values v1 and v2 for STL formulae phi1 and
% phi2, compute alpha(v1, v2) according to the name of the semantics 
% requested.

% alpha is sound when, for x1 and x2 >= 0, alpha(x1, x2) > 0 only if 
% x1 > 0 and x2 > 0. 0 otherwise.

function res = s_alpha(v1, v2, semantics)
    switch semantics
        case 'max-breach'
            res = min(v1, v2);
        case 'const-breach'
            res = min(v1, v2);
        case 'plus-breach'
            res = koenAndPlus(v1, v2);
        case 'telex'
            res = min(v1, v2);
        case 'belta'
            res = min(v1, v2);
        case 'agm-product'
            res = safemult((1+v1),(1+v2))-1;
        case 'sum-product'
            res = safemult(v1, v2);
        case 'sum-min'
            res = min(v1, v2);
        case 'max-product'
            res = safemult(v1, v2);
        case 'minonly'
            res = min(v1, v2);
        case 'smoothrect'
            res = min(v1, v2);
        case 'smooth1'
            res = safemult(v1, v2);
        case 'smooth2'
            %res = log(exp(v1) + exp(v2));
            res = v1 + v2;
        otherwise
            error('Unknown semantics for alpha!');
    end
end