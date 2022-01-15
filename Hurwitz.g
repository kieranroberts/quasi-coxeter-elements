swap := function(set, first, second)
    local temp;

    temp := set[first];
    set[first] := set[second];
    set[second] := temp; 
    
    return true;
end;


# f is a function of i that returns the braid generator sigma_i 
# as a function.

f := function(i)
    return function(g)
        local gCopy;

        gCopy := MutableCopyMat(g);
        swap(gCopy,i,i+1);
        gCopy[i+1] := gCopy[i+1]^gCopy[i];
       return gCopy;
    end;
end;



# invf is a function of i that returns the braid generator sigma_i^{-1}
# as a function of a tuple g.

invf := function(i)
    return function(g)
        local gCopy;

        gCopy := MutableCopyMat(g);
        swap(gCopy,i,i+1);
        gCopy[i] := gCopy[i]^gCopy[i+1];
       return gCopy;
    end;
end;


#########################################################################
# SigmaAction constructs the braid generators and the inverses and
# places them into lists sigma and inv_sigma, respectively.
#
# sigma = [sigma[1], sigma[2], ... , sigma[k-1]]
# inv_sigma = [ sigma[1]^-1, sigma[2]^-1, ... , sigma[k-1]^-1]
#########################################################################
 
SigmaAction := function(sigma, inv_sigma,k) 
    local i, sig, inv_sig;

    sig:=[]; inv_sig := [];

    # computes the sigma function for each i
    for i in [1..k-1] do
        sig[i] := f(i);
        Add(sigma, sig[i]);
    od;
   
    # computes the inverse sigma function for each i
    for i in [1..k-1] do
        inv_sig[i] := invf(i);
        Add(inv_sigma, inv_sig[i]);
    od;
    return true;
end;


#############################################################################
# BraidAction defines the braid action on groups.
#
# 1. For a given group element g = t_1 ... t_n. By an abuse of notation 
#    we write g = [t_1, ... , t_n]. 
# 2. We take an element s from the braid group B(n) and construct 
#    the braid generators [sigma[1], ... , sigma[n-1]] and its inverses.
# 3. We break down s = s_1 s_2 ... s_k, where each s_i is a braid generator
#    or its inverse and apply them to the tuple g = [t_1, ..., t_n] from 
#    right to left.
#
#############################################################################

BraidAction:=function(g,s)  
    local i,j,k, sigma, inv_sigma, l, e;
    
    sigma := []; inv_sigma :=[];
    SigmaAction(sigma, inv_sigma, Length(g));
    l := NumberSyllables(s);
    
    for i in [1..l] do
        j:=l+1-i;
        e := ExponentSyllable(s,j); 
        if e>0 then
            for k in [1..e] do
                g:=sigma[GeneratorSyllable(s,j)](g);
            od;
            else
                for k in [1..AbsInt(e)] do
                    g:=inv_sigma[GeneratorSyllable(s,j)](g);
                od;  
       fi;
    od;
    return g;
end;

###############################################################################
#
#  Hurwitz Orbit computes the orbit of a tuple of reflections under the 
#  braid action.
#
###############################################################################

HurwitzOrbit := function(wt)
    local j;
    
    j := Size(wt);
    return Orbit(FreeGroup(j-1), wt, BraidAction);
end;
