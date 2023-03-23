function res = beltamax(v1,v2)
    %arbitrarily choosing beta=10 
    beta=10;
    res = (1/beta)*log(exp(beta*v1)+exp(beta*v2));
end