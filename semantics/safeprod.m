function res = safeprod(v)
% Takes the product of a vector, but treats NaN's as zeros.
    res = prod(v);
    res(isnan(res)) = 0;
end