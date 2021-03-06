##################################################################################################
# This file will contain generating sets for each of the irreducible Coxeter groups A_n (n>=1), 
# B_n (n>=2), D_n (n>=4), E_6, E_7, E_8, F_4, G_2, H_3, H_4. The generators along with the groups
# will give the usual Coxeter systems. Furthermore, the permutation representations of the groups 
# are minimal.
#
##################################################################################################


##################################################################################################
# Infinite Families A_n, B_n, D_n:
# W(A_n) = Sym(n+1) with fundamental reflections { (1,2), (2,3), ... , (n, n+1) }.
# W(B_n) = Sym(\pm n) is a permutation group on 2n letters.
# W(D_n) is an index-2 subgroup of W(B_n) that consists of the "even" permutations of W(B_n).
#
##################################################################################################

CoxGrpGensA := function(n)
    local i;
    return List([1..n],i->(i,i+1));
end;

CoxGrpGensB := function(n)
    local i;
    return Union(List([1..n-1], i->(i,i+1)(i+n,(i+1)+n)),[(1,n+1)]);
end;

CoxGrpGensD := function(n)
    local i;
    return Union(List([1..n-1], i->(i,i+1)(i+n,(i+1)+n)),[(1,n+2)(2,n+1)]);
end;



##################################################################################################
# We adopt the Bourbaki diagrams for the exceptional Coxeter groups E_6, E_7, E_8:
#
# 1 -- 3 -- 4 -- 5 -- 6
#           |
#           2
#
# 1 -- 3 -- 4 -- 5 -- 6 -- 7
#           |
#           2
#
# 1 -- 3 -- 4 -- 5 -- 6 -- 7 -- 8
#           |
#           2
#
##################################################################################################

CoxGrpGensE := function(n)
    if n = 6 then
        return [ \
                (1,2)(8,12)(9,13)(10,14)(11,15)(22,23), \
                (1,12)(2,8)(3,7)(19,27)(20,26)(21,25), \
                (2,3)(7,8)(13,16)(14,17)(15,18)(23,24), \
                (3,4)(8,9)(12,13)(17,19)(18,20)(24,25), \
                (4,5)(9,10)(13,14)(16,17)(20,21)(25,26), \
                (5,6)(10,11)(14,15)(17,18)(19,20)(26,27) \
              ];
    elif n = 7 then
        return [ \
                (1,2)(5,14)(7,21)(10,18)(15,20)(17,28)(29,30), \
                (1,18)(2,10)(3,9)(4,8)(13,25)(24,26)(29,30), \
                (2,9)(3,10)(5,23)(11,28)(15,16)(21,27)(29,30), \
                (1,15)(2,20)(3,22)(8,27)(11,13)(23,26)(29,30), \
                (4,13)(5,15)(8,25)(12,22)(14,20)(16,23)(29,30), \
                (4,24)(5,21)(6,12)(7,14)(8,26)(23,27)(29,30), \
                (4,25)(6,19)(7,17)(8,13)(11,27)(21,28)(29,30) \
              ];
    elif n = 8 then
        return [ \
                   (1,121)(3,9)(11,16)(17,23)(19,24)(25,30)(27,31)(32,37)(33,38)(35,39)(40,45)(41,46)(43,47)(48,52)(49,53) \
                   (50, 54)(55, 59)(56,60)(61, 66)(62, 67)(68,73)(74, 79)(93, 97)(98, 101)(102, 104)(105, 107)(108, 109)(110,111)(112, 113) \
                   (123, 129)(131, 136)(137, 143)(139, 144)(145, 150)(147,151)(152, 157)(153,158)(155, 159)(160, 165)(161, 166)(163, 167) \
                   (168, 172)(169,173)(170, 174)(175, 179)(176, 180)(181, 186)(182, 187)(188, 193)(194,199)(213, 217)(218, 221)(222,224) \
                   (225, 227)(228, 229)(230, 231)(232, 233),\
                   \
                    (2, 122)(4, 10)(11, 17)(12, 18)(16, 23)(19, 25)(20, 26)(24, 30)(27,33)(28, 34)(31, 38)(35, 41)(36, 42)(39, 46)(43, 50)\
                   (47, 54)(63, 69)(70,75)(76, 80)(77, 81)(82,85)(83, 86)(87, 90)(88, 91)(92, 95)(96, 100)(110, 112)(111,113)(114, 115)\
                   (124, 130)(131, 137)(132, 138)(136, 143)(139, 145)(140,146)(144, 150)(147, 153)(148,154)(151, 158)(155, 161)(156, 162)\
                   (159, 166)(163, 170)(167,174)(183, 189)(190, 195)(196, 200)(197, 201)(202, 205)(203, 206)(207,210)(208, 211)(212, 215)\
                   (216,220)(230, 232)(231, 233)(234, 235),\
                   \
                  (1, 9)(3, 123)(4, 11)(10, 17)(12, 19)(18, 25)(20, 27)(26, 33)(28,35)(34, 41)(36, 43)(37, 44)(42, 50)(45, 51)(52, 57)\
                  (53, 58)(59, 64)(60,65)(66, 71)(67, 72)(73,78)(79, 84)(89, 93)(94, 98)(99, 102)(103, 105)(106, 108)(111,114)(113, 115)\
                  (121, 129)(124, 131)(130, 137)(132, 139)(138, 145)(140,147)(146, 153)(148,155)(154, 161)(156, 163)(157, 164)(162, 170)\
                  (165, 171)(172,177)(173, 178)(179, 184)(180, 185)(186,191)(187, 192)(193, 198)(199,204)(209, 213)(214, 218)(219,222)\
                  (223, 225)(226, 228)(231, 234)(233, 235),\
                  \
                  (2, 10)(3, 11)(4, 124)(5, 12)(9, 16)(13, 20)(21, 28)(25, 32)(29,36)(30, 37)(33, 40)(38, 45)(41, 49)(46, 53)(50, 56)\
                  (54, 60)(57, 63)(64,70)(71, 76)(72, 77)(78,83)(84, 88)(85, 89)(90, 94)(95, 99)(100,103)(108, 110)(109,111)(115, 116)\
                  (122, 130)(123, 131)(125, 132)(129, 136)(133, 140)(141,148)(145, 152)(149, 156)(150,157)(153, 160)(158, 165)(161, 169)\
                  (166, 173)(170, 176)(174,180)(177, 183)(184, 190)(191, 196)(192, 197)(198, 203)(204, 208)(205,209)(210, 214)(215, 219)\
                  (220,223)(228, 230)(229, 231)(235, 236),\
                  \
                 (4, 12)(5, 125)(6, 13)(10, 18)(11, 19)(14, 21)(16, 24)(17, 25)(22,29)(23, 30)(40, 48)(45, 52)(49, 55)(51, 57)(53, 59)\
                 (56, 62)(58, 64)(60,67)(65, 72)(76, 82)(80,85)(83, 87)(86, 90)(88, 92)(91, 95)(103, 106)(105, 108)(107,109)(116, 117)\
                 (124, 132)(126, 133)(130, 138)(131, 139)(134, 141)(136,144)(137, 145)(142, 149)(143,150)(160, 168)(165, 172)(169, 175)\
                 (171, 177)(173, 179)(176,182)(178, 184)(180, 187)(185, 192)(196, 202)(200, 205)(203, 207)(206,210)(208, 212)(211, 215)\
                 (223,226)(225, 228)(227, 229)(236, 237),\
                 \
                (5, 13)(6, 126)(7, 14)(12, 20)(15, 22)(18, 26)(19, 27)(24, 31)(25,33)(30, 38)(32, 40)(37, 45)(44, 51)(55, 61)(59, 66)(62, 68)\
                (64, 71)(67,73)(70, 76)(72, 78)(75,80)(77, 83)(81, 86)(92, 96)(95, 100)(99,103)(102, 105)(104,107)(117, 118)(125, 133)\
                (127, 134)(132, 140)(135, 142)(138, 146)(139,147)(144, 151)(145, 153)(150,158)(152, 160)(157, 165)(164, 171)(175, 181)\
                (179, 186)(182,188)(184, 191)(187, 193)(190, 196)(192, 198)(195, 200)(197, 203)(201,206)(212, 216)(215, 220)(219,223)(222, 225)\
                (224, 227)(237, 238),\
                \
               (6, 14)(7, 127)(8, 15)(13, 21)(20, 28)(26, 34)(27, 35)(31, 39)(33,41)(38, 46)(40, 49)(45, 53)(48, 55)(51, 58)(52, 59)(57, 64)\
               (63, 70)(68,74)(69, 75)(73, 79)(78,84)(83, 88)(86, 91)(87, 92)(90, 95)(94,99)(98, 102)(101,104)(118, 119)(126, 134)(128, 135)\
               (133, 141)(140, 148)(146, 154)(147,155)(151, 159)(153, 161)(158,166)(160, 169)(165, 173)(168, 175)(171, 178)(172, 179)(177,184)\
               (183, 190)(188, 194)(189, 195)(193, 199)(198, 204)(203, 208)(206,211)(207, 212)(210, 215)(214,219)(218, 222)(221, 224)(238, 239),\
               (7, 15)(8, 128)(14, 22)(21, 29)(28, 36)(34, 42)(35, 43)(39, 47)(41,50)(46, 54)(49, 56)(53, 60)(55, 62)(58, 65)(59, 67)(61, 68)\
               (64, 72)(66,73)(70, 77)(71, 78)(75,81)(76, 83)(80, 86)(82, 87)(85, 90)(89, 94)(93, 98)(97,101)(119, 120)(127, 135)(134, 142)\
               (141, 149)(148, 156)(154, 162)(155,163)(159, 167)(161, 170)(166,174)(169, 176)(173, 180)(175, 182)(178, 185)(179, 187)(181,188)\
               (184, 192)(186, 193)(190, 197)(191, 198)(195, 201)(196, 203)(200,206)(202, 207)(205, 210)(209,214)(213, 218)(217, 221)(239, 240)];
    else
        Error("<n> must be 6,7,8");
    fi;
end;

##################################################################################################
# The remaining groups are F_4, G_2, H_3, H_4
#
##################################################################################################

CoxGrpGensF := function()
    return [ (1,2)(4,5)(7,8)(10,11)(13,19)(14,21)(16,22)(20,24)(18,23), \
            (2,3)(5,6)(8,9)(11,12)(13,14)(15,16)(17,18)(19,20)(21,24), \
            (3,6)(9,12)(14,16)(18,20)(21,22)(23,24), \
            (1,22)(2,16)(3,15)(10,23)(11,18)(12,17)];
end;

CoxGrpGensG := function()
    return [ (1,3)(4,6), (1,4)(2,3)(5,6) ];
end;

CoxGrpGensH := function(n)
    if n = 3 then
        return [ (1,3)(4,5)(6,7), \
                (1,5)(2,4)(6,7), \
                (1,5)(3,4)(6,7) ];
    elif n = 4 then
        return [(1,2)(3,4)(5,7)(8,11)(10,14)(12,16)(15,19)(17,22)(20,26)(21,27)(23,30)(28,34)(29,35)(33,41)(36,43)(37,44)(38,47)(40,49)(42,51)\
                (45,53)(46,56)(48,58)(52,61)(54,64)(55,65)(59,69)(62,72)(63,74)(66,77)(68,80)(70,82)(73,85)(75,87)(76,89)(78,91)(81,93)(83,96)\
                (86,97)(88,98)(90,99)(92,102)(94,104)(95,105)(100,108)(101,109)(106,112)(107,114)(110,115)(113,117)(119,120),\
                \
                (2,3)(4,6)(7,10)(9,12)(11,15)(13,17)(14,18)(16,21)(19,25)(22,29)(24,28)(26,33)(30,38)(31,36)(32,40)(34,42)(37,46)(39,48)(41,50)\
                (43,52)(44,54)(45,55)(47,57)(49,59)(53,63)(56,66)(58,68)(60,70)(62,73)(64,75)(65,76)(67,78)(71,83)(72,84)(74,86)(77,87)(79,92) \
                (81,90)(82,95)(89,97)(91,101)(93,103)(94,100)(96,106)(98,107)(102,110)(104,111)(113,116)(114,118)(117,119),\
                \
                (3,5)(4,7)(6,9)(10,12)(14,16)(15,20)(17,23)(18,24)(19,26)(21,28)(22,30)(25,32)(27,34)(29,37)(31,39)(33,40)(35,44)(36,45)(38,46)\
                (41,49)(43,53)(47,56)(48,55)(50,60)(52,62)(57,67)(58,65)(59,70)(61,72)(63,73)(66,78)(68,81)(69,82)(71,79)(74,85)(75,88)(76,90)\
                (77,91)(80,93)(83,94)(87,98)(89,99)(92,100)(96,104)(101,107)(102,108)(106,113)(109,114)(111,116)(112,117),\
                \
                (5,8)(7,11)(9,13)(10,15)(12,17)(14,19)(16,22)(18,25)(20,23)(21,29)(24,31)(26,30)(27,35)(28,36)(32,39)(33,38)(34,43)(37,45)(40,48)\
                (41,47)(42,52)(44,53)(46,55)(49,58)(50,57)(51,61)(54,63)(56,65)(59,68)(60,71)(64,74)(66,76)(67,79)(69,80)(70,83)(75,86)(77,89)(78,92)\
                (81,94)(82,96)(87,97)(90,100)(91,102)(93,104)(95,106)(99,108)(101,110)(103,111)(105,112)(109,115)];
    else
        Error("<n> must be 3,4");
    fi;
end;
