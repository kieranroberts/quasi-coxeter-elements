###############################################################################
#
# QuasiCoxeterClasses(W,T,n,q)
#
# W is the Coxeter group
# T is the set of reflections
# n is the rank of the Coxeter group
# q is the number of conjugacy classes of quasi-Coxeter elements which you can 
#   look up in noqcclasses.txt
#
# The return value is a list of pairs (w, x) where w is a representative from
# the conjugacy class and x is a set of n reflections whose product gives w.
#
# So for example if you define W to be D_4 and T to be its reflections. Then 
# typing:
# qcc := QuasiCoxeterClasses(W,T,4,2);
# will produce a list of of two pairs [ [w1, wt1], [w2, wt2] ] and one of w1
# or w2 will be the Coxeter element and the other will be the proper quasi-
# Coxeter element.
#
############################################################################## 


RandomNTuple := function(T,n)
    local tuple, i, x;
    
    tuple := [];
    
    for i in [ 1 .. n ] do
        x := Random(T);
        while x in tuple do
            x := Random(T);
        od;
        Add(tuple, x);
    od;
    
    return tuple;
end;

IsQCE := function(W, redw)
    return (W = Group(redw));
end;

IsConj := function(W, S, w)
    local flag, s;
    
    flag := false;
    for s in S do
        if IsConjugate(W, s, w) then
            flag := true;
        fi;
   od;
   return flag;
end;


QuasiCoxeterClasses := function(W,T,n,q)
    local QCox, QC, t, wt;

    QCox := [];
    QC := [];
    
    while Size(QCox) < q do
        t := RandomNTuple(T,n);
        wt := Product(t);
        while IsConj(W, QCox, wt) or not(IsQCE(W,t)) do
            t := RandomNTuple(T,n);
            wt := Product(t);
        od;
        Add(QCox, wt);
        Add(QC,[wt,t]);
    od;
    return QC;
end;
