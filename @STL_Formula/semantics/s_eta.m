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