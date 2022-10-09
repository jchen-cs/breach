% The nu function, which implements the rectifier for rho-
% For a vector of robustness scores v and times t, compute the result of 
% mu(v) according to the name of the semantics requested. 

% mu is sound when mu(x) < 0 only if x < 0, and 0 otherwise.

function res = s_mu(t, v, semantics)
    switch semantics
        case 'max-breach'
            res = min(v, 0);
        case 'const-breach'
            res = min(-100*sign(v), 0);
        case 'plus-breach'
            res = min(v, 0);
        case 'telex'
            res = min(TeLExPeak(v, 1), 0);
        case 'belta'
            res = min(v, 0);
        case 'agm-product'
            res = min(v, 0);
        case 'sum-product'
            res = min(v, 0);
        case 'sum-min'
            res = min(v, 0);
        case 'max-product'
            res = min(v, 0);
        case 'minonly'
            res = min(v, 0);
        case 'smoothrect'
            res = v;
            res(res >= 0) = 0;
            res(res < 0) = res(res < 0) .* exp(1 ./ res(res < 0));
        otherwise
            error('Unknown semantics for mu!');
    end
end