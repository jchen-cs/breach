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
