function res = TeLExExpand(gamma, t1, t2)
    res = 2 / (1+exp(-gamma*(t2-t1+1)));
end