function res = safemult(v1, v2)
% Takes the element-wise product of two vectors, but treats NaN's as zeros.
    res = v1 .* v2;
    res(isnan(res)) = 0;
end