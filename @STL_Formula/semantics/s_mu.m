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
            res(res < 0) = res .* exp(1 ./ res);
        otherwise
            error('Unknown semantics for mu!');
    end
end