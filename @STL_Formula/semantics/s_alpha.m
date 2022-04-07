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
            error('Unknown semantics for alpha!');
    end
end