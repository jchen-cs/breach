function res = beltamaxvector(v)
    %arbitrarily choosing beta=10 
    beta=10;
    res = (1/beta)*log(sum(exp(beta*v)));
end