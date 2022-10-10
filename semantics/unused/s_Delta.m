% The Delta function, which is the pessimistic trace integrator typically 
% used for UNTIL in rho+

% That is, it takes the "worst" score within a trace.

% For a vector of robustness values v, compute Delta(v) according to the 
% name of the semantics requested.

% Delta is sound when, for robustness values x where all x >= 0, 
% Delta(x) > 0 only if all x > 0.

function res = s_Delta(v, semantics)
    switch semantics
        case 'max-breach'
            res = min(v);
        case 'const-breach'
            res = min(v);
        case 'plus-breach'
            res = min(v);
        case 'telex'
            res = min(v);
        case 'belta'
            res = min(v);
        case 'agm-product'
            res = prod(v+1)-1;
        case 'sum-product'
            res = prod(v);
        case 'sum-min'
            res = min(v);
        case 'max-product'
            res = prod(v);
        case 'minonly'
            res = min(v);
        case 'smoothrect'
            res = min(v);
        otherwise
            error('Unknown semantics for Delta!');
    end
end
