Read("Hurwitz.g");

###############################################################################
#
#  Utility functions refKProduct(s,k) and refKtuple(s,k) that returns all 
#  k-tuples and products of k-tuples formed from the k-prefixes of the set s 
#  of tuples  
#
###############################################################################

refKProduct := function(s,k)
    local x, prd;
    
    prd := [];
    for x in s do
        Add(prd, Product(x{[1..k]}));
    od;
    return Set(prd);
end;

refKTuple := function(s,k)
    local x, tup;
    
    tup := [];
    for x in s do
        Add(tup, x{[1..k]});
    od;
    return Set(tup);
 end;  


###############################################################################
#
#  ncp := NCP(wt) takes a tuple of reflections, which are assumed to be a 
#  reduced expression of a quasi-Coxeter element and returns the noncrossing 
#  partition of w = Product(wt) and returns a list with 5 elements 
#  corresponding to the 5 layers:
#
#  ncp[i] is layer i-1 and consists of the pairs:
#  [ g, [ [ t_1 ... t_{i-1} ] : g = Product(t_1 ... t_{i-1} ) ] ]
#
###############################################################################

NCP := function(wt)
    local s,n, i, prods, tups, res, x, temp1, temp2, g;
    
    s := HurwitzOrbit(wt);
    n := Length(wt);
    
    res := [[ [(), [[()]] ]]];

    for i in [1..n] do
        prods := refKProduct(s,i);
        tups := refKTuple(s,i);
        temp1 := [];
        for g in prods do
            temp2 := [];
            for x in tups do
                if g = Product(x) then
                    Add(temp2, x);
                fi;
            od;
            Add(temp1, [g,temp2]);
        od;
        Add(res, temp1);
    od;
    return res;
end;
