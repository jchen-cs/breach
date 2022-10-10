%% Matches two traces in time, using interpolation/extrapolation to fill in missing samples
function [t, v1, v2] = matchTraces(time_values1, time_values2, valarray1, valarray2)
    t = union(time_values1, time_values2);
    v1 = interp1(time_values1, valarray1, t, 'previous', 'extrap');
    v2 = interp1(time_values2, valarray2, t, 'previous', 'extrap');
end