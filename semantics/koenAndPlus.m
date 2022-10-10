function res = koenAndPlus(v1, v2)
    % Applies Koen's &+ integrator
    % Default is just min
    res = min(v1, v2);
    % If both are negative, add them.
    both_negative = v1 < 0 & v2 < 0;
    res(both_negative) = v1(both_negative) + v2(both_negative);
    % If both are positive, perform "parallel resistance" addition
    both_positive = v1 > 0 & v2 > 0;
    res(both_positive) = 1 ./ ((1./v1(both_positive)) + (1./v2(both_positive)));
end