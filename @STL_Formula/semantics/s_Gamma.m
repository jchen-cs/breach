function res = s_Gamma(v, semantics, interval)
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
            res = sum(v);
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
            error('Unknown semantics for Gamma!');
    end
end