% The zeta function, which is the pessimistic integrator typically used
% for UNTIL in rho+
% That is, it takes the "worse" score of its two inputs.

% For vectors of robustness values v1 and v2 for STL formulae phi1 and
% phi2, compute zeta(v1, v2) according to the name of the semantics 
% requested.

% zeta has the same soundness conditions as alpha.

function res = s_zeta(v1, v2, semantics)
    switch semantics
        case 'max-breach'
            res = min(v1, v2);
        case 'const-breach'
            res = min(v1, v2);
        case 'plus-breach'
            res = min(v1, v2);
        case 'telex'
            res = min(v1, v2);
        case 'belta'
            res = min(v1, v2);
        case 'agm-product'
            res = ((1+v1).*(1+v2))-1;
        case 'sum-product'
            res = v1 .* v2;
        case 'sum-min'
            res = min(v1, v2);
        case 'max-product'
            res = v1 .* v2;
        case 'minonly'
            res = min(v1, v2);
        case 'smoothrect'
            res = min(v1, v2);
        otherwise
            error('Unknown semantics for zeta!');
    end
end