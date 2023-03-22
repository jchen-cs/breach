function res = smooth2alpha(v1,v2)
    % Get a logical array as long as both input vectors, with a 1 in any
    % index where both v1 and v2 are positive
    isNonZero = (v1 > 0 & v2 > 0);
    % Initialize the results vector with all zeros (the default answer)
    res = zeros(size(isNonZero));
    res(isNonZero) = (1./(log(0.5*(exp(1./v1(isNonZero))+exp(1./v2(isNonZero))))));
end

