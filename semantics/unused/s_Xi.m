% The Xi function, which is the optimistic trace integrator typically 
% used for UNTIL in rho-

% That is, it takes the "best" score within a trace.

% For a vector of robustness values v, compute Gamma(v) according to the 
% name of the semantics requested.

% Xi is sound when, for robustness values x where all x >= 0, 
% Xi(x) > 0 only if there exists at least one x such that x > 0.

% When used in rho-, Xi is usually negated (ie -Xi(-v1, -v2) ) so 
% it's treated as an pessimistic integrator.

function res = s_Xi(v, semantics, interval)
    switch semantics
        case 'max-breach'
            res = max(v);
        case 'const-breach'
            res = max(v);
        case 'plus-breach'
            res = max(v);
        case 'telex'
            res = TeLEXExpand(0.01, interval(1), interval(2)) .* max(v);
        case 'belta'
            res = max(v);
        case 'agm-product'
            res = max(v);
        case 'sum-product'
            res = sum(v);
        case 'sum-min'
            res = sum(v);
        case 'max-product'
            res = max(v);
        case 'minonly'
            res = min(v);
        case 'smoothrect'
            res = max(v);
        otherwise
            error('Unknown semantics for Xi!');
    end
end