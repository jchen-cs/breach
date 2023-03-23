function [time_values, valarray_P, valarray_N] = PlusMinusUntil(time_values1, valarray_P1, valarray_N1, time_values2, valarray_P2, valarray_N2, interval, I___, semantics) 
    switch semantics
        case 'max-breach'
            l_zeta = @min;
        case 'const-breach'
            l_zeta = @min;
        case 'plus-breach'
            l_zeta = @min;
        case 'telex'
            l_zeta = @min;
        case 'belta'
            %l_zeta = @min;
            l_zeta = @(v1, v2)(-beltamax(-v1, -v2));
        case 'agm-product'
            l_zeta = @(v1, v2)(safemult((1+v1),(1+v2)))-1;
        case 'sum-product'
            l_zeta = @safemult;
        case 'sum-min'
            l_zeta = @min;
        case 'max-product'
            l_zeta = @safemult;
        case 'minonly'
            l_zeta = @min;
        case 'smoothrect'
            l_zeta = @min;
        case 'smooth1'
            l_zeta = @safemult;
        case 'smooth2'
            % This seems to cause problems
            %l_zeta = @(v1, v2)(-log(exp(-v1) + exp(-v2)));
            l_zeta = @smooth2alpha;
        otherwise
            error('Unknown semantics for zeta!');
    end
    
    switch semantics
        case 'max-breach'
            l_eta = @max;
        case 'const-breach'
            l_eta = @max;
        case 'plus-breach'
            l_eta = @max;
        case 'telex'
            l_eta = @max;
        case 'belta'
            %l_eta = @max;
            l_eta = @beltamax;
        case 'agm-product'
            l_eta = @max;
        case 'sum-product'
            l_eta = @(v1, v2)(v1 + v2);
        case 'sum-min'
            l_eta = @(v1, v2)(v1 + v2);
        case 'max-product'
            l_eta = @max;
        case 'minonly'
            l_eta = @min;
        case 'smoothrect'
            l_eta = @max;
        case 'smooth1'
            l_eta = @(v1, v2)(v1 + v2);
        case 'smooth2'
            l_eta = @(v1, v2)log((exp(v1) + exp(v2))/2);
        otherwise
            error('Unknown semantics for eta!');
    end
    
    switch semantics
        case 'max-breach'
            l_Gamma = @(v, interval)max(v);
        case 'const-breach'
            l_Gamma = @(v, interval)max(v);
        case 'plus-breach'
            l_Gamma = @(v, interval)max(v);
        case 'telex'
            l_Gamma = @(v, interval)(safemult(TeLExExpand(0.01, interval(1), interval(2)), max(v)));
        case 'belta'
            l_Gamma = @(v, interval)sum(v);
        case 'agm-product'
            l_Gamma = @(v, interval)max(v);
        case 'sum-product'
            l_Gamma = @(v, interval)sum(v);
        case 'sum-min'
            l_Gamma = @(v, interval)sum(v);
        case 'max-product'
            l_Gamma = @(v, interval)max(v);
        case 'minonly'
            l_Gamma = @(v, interval)min(v);
        case 'smoothrect'
            l_Gamma = @(v, interval)max(v);
        case 'smooth1'
            l_Gamma = @(v, interval)sum(v);
        case 'smooth2'
            l_Gamma = @(v, interval)(log(mean(exp(v))));
            %l_Gamma = @(v, interval)sum(v);
        otherwise
            error('Unknown semantics for Gamma!');
    end
    
    switch semantics
        case 'max-breach'
            l_Delta = @min;
        case 'const-breach'
            l_Delta = @min;
        case 'plus-breach'
            l_Delta = @min;
        case 'telex'
            l_Delta = @min;
        case 'belta'
            %l_Delta = @min;
            l_Delta = @(v)(-beltamaxvector(-v));
        case 'agm-product'
            l_Delta = @(v)(safeprod(v+1)-1);
        case 'sum-product'
            l_Delta = @safeprod;
        case 'sum-min'
            l_Delta = @min;
        case 'max-product'
            l_Delta = @safeprod;
        case 'minonly'
            l_Delta = @min;
        case 'smoothrect'
            l_Delta = @min;
        case 'smooth1'
            l_Delta = @safeprod;
        case 'smooth2'
            %l_Delta = @(v)(-log(sum(exp(-v))));
            l_Delta = @(v)(1/log(mean(exp(1./v))));
        otherwise
            error('Unknown semantics for Delta!');
    end
    
    switch semantics
        case 'max-breach'
            l_Xi = @(v, interval)max(v);
        case 'const-breach'
            l_Xi = @(v, interval)max(v);
        case 'plus-breach'
            l_Xi = @(v, interval)max(v);
        case 'telex'
            l_Xi = @(v, interval)(safemult(TeLExExpand(0.01, interval(1), interval(2)), max(v)));
        case 'belta'
            %l_Xi = @(v, interval)max(v);
            l_Xi = @(v, interval)beltamaxvector(v);
        case 'agm-product'
            l_Xi = @(v, interval)max(v);
        case 'sum-product'
            l_Xi = @(v, interval)sum(v);
        case 'sum-min'
            l_Xi = @(v, interval)sum(v);
        case 'max-product'
            l_Xi = @(v, interval)max(v);
        case 'minonly'
            l_Xi = @(v, interval)min(v);
        case 'smoothrect'
            l_Xi = @(v, interval)max(v);
        case 'smooth1'
            l_Xi = @(v, interval)sum(v);
        case 'smooth2'
            l_Xi = @(v, interval)(log(mean(exp(v))));
        otherwise
            error('Unknown semantics for Xi!');
    end
    
    switch semantics
        case 'max-breach'
            l_Theta = @min;
        case 'const-breach'
            l_Theta = @min;
        case 'plus-breach'
            l_Theta = @min;
        case 'telex'
            l_Theta = @min;
        case 'belta'
            %l_Theta = @min;
            l_Theta = @sum;
        case 'agm-product'
            l_Theta = @(v)(safeprod(v+1)-1);
        case 'sum-product'
            l_Theta = @safeprod;
        case 'sum-min'
            l_Theta = @min;
        case 'max-product'
            l_Theta = @safeprod;
        case 'minonly'
            l_Theta = @min;
        case 'smoothrect'
            l_Theta = @min;
        case 'smooth1'
            l_Theta = @safeprod;
        case 'smooth2'
            %l_Theta = @(v)(-log(sum(exp(-v))));
            l_Theta = @(v)(1/log(mean(exp(1./v))));
        otherwise
            error('Unknown semantics for Theta!');
    end
    
    % If time horizons are finite, duplicate the final values at the end.
    if(I___(end)~=inf)
        time_values1 = [time_values1 time_values1(end)+I___(end)];
        valarray_P1 = [valarray_P1 valarray_P1(end)];
        valarray_N1 = [valarray_N1 valarray_N1(end)];
        time_values2 = [time_values2 time_values2(end)+I___(end)];
        valarray_P2 = [valarray_P2 valarray_P2(end)];
        valarray_N2 = [valarray_N2 valarray_N2(end)];
    end
    
    % Length-match the traces
    [time_values_combined, P1, P2] = matchTraces(time_values1, time_values2, valarray_P1, valarray_P2);
    [~, N1, N2] = matchTraces(time_values1, time_values2, valarray_N1, valarray_N2);
    
    % Get times of all samples between 1 and 2
    
    % We're interested in time values inside the interval of interest
    time_values = time_values_combined((time_values_combined >= interval(1)) & (time_values_combined <= interval(2)));
    valarray_P = zeros(size(time_values));
    valarray_N = zeros(size(time_values));
    N = size(time_values, 2);
    assert(all(valarray_P1 >= 0));
    assert(all(valarray_P2 >= 0));
    assert(all(valarray_N1 <= 0));
    assert(all(valarray_N2 <= 0));
    assert(all(~((valarray_P1 > 0) & (valarray_N1 < 0))))
    assert(all(~((valarray_P2 > 0) & (valarray_N2 < 0))))
    parfor k = 1:N
        current_time = time_values(k);
        time_start = current_time + I___(1);
        time_end = current_time + I___(2);
        idx_start = find(time_values_combined >= time_start, 1);
        idx_end = find(time_values_combined >= time_end, 1);
        % For each time inside the Until interval starting from time at k: 
        zeta_result = zeros(1, idx_end - idx_start);
        eta_result = zeros(1, idx_end - idx_start);
        % k_plus_kprime represents k+k' in the semantics, or the indices of
        % time within the Until interval (which is a positive offset from
        % the current time index, k)
        % Here we look at all possible times for phi to "activate" within
        % the interval, and find the worst-case.
        
        for k_plus_kprime = idx_start:idx_end
            % For psi Until phi (here, trace 1 Until trace 2), look at 
            % psi's robustness from current time until the time phi takes 
            % effect (k_plus_kprime)
            psi_indices = k:k_plus_kprime;
            % Delta integrates all of trace psi between t(k) and t(k_prime)
            % Same for Xi but for negative robustness
            
            % zeta integrates the result of Delta, as well as the
            % robustness value of trace 2 at k+k_prime
            % This gives us lists of values, where the i-th value in each
            % list corresponds to zeta/eta for each value of k'
            zeta_result(k_plus_kprime - idx_start + 1) = l_zeta(P2(k_plus_kprime), l_Delta(P1(psi_indices)));
            eta_result(k_plus_kprime - idx_start + 1) = l_eta(-N2(k_plus_kprime), l_Xi(-N1(psi_indices), I___));
        end
        valarray_P(k) = l_Gamma(zeta_result, I___);
        valarray_N(k) = -l_Theta(eta_result); % Minus?

        % Soundness check: Theta and Delta have the same requirements:
        % If all x_k >= 0 and Delta(x_k) > 0 then \forall k. x_k > 0
        % DANGER! Disabling soundness assertions can cause unexpected
        % behavior!
        % assert(valarray_N(k) <= 0)
        % assert(valarray_P(k) >= 0)
        
    end
    % assert(all(valarray_N <= 0))
    % assert(all(valarray_P >= 0))
    % Assert that valarray_P and valarray_N are mutually exclusive, ie you
    % can't have both be nonzero
    % assert(all(~((valarray_P ~= 0) & (valarray_N ~= 0))))
    %size(valarray_P)
end
