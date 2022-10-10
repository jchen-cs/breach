function res = TeLExPeak(r, beta)
    res = (1 ./(r+exp(-beta .* r))) - exp(-r);
end