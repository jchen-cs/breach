% The nu function, which implements the rectifier for rho+
% For a vector of robustness scores v and times t, compute the result of 
% nu(v) according to the name of the semantics requested. 

% nu is sound when nu(x) > 0 only if x > 0, and 0 otherwise.

function res = s_nu(t, v, semantics)
    switch semantics
        case 'max-breach'
            res = max(v, 0);
        case 'const-breach'
            res = max(100*sign(v), 0);
        case 'plus-breach'
            res = max(v, 0);
        case 'telex'
            res = max(TeLExPeak(v, 1), 0);
        case 'belta'
            res = max(v, 0);
        case 'agm-product'
            res = max(v, 0);
        case 'sum-product'
            res = max(v, 0);
        case 'sum-min'
            res = max(v, 0);
        case 'max-product'
            res = max(v, 0);
        case 'minonly'
            res = max(v, 0);
        case 'smoothrect'
            res = v;
            res(res <= 0) = 0;
            res(res > 0) = res(res > 0) .* exp(-1 ./ res(res > 0));
        case 'smooth1'
            res = max(v .* exp(-1 ./ v), 0);
        case 'smooth2'
            res = -max(-v .* exp(1 ./ v), 0);
        otherwise
            error('Unknown semantics for nu!');
    end
end