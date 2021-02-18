# fix for ice system connectivity issue

# Fix H
#  I=     6 IAn=  1 IBond, RBType=   326  1.00   347  0.10
#  I=     7 IAn=  1 IBond, RBType=   327  1.00
#                                            ^
# sed -r -i '/IAn=  1/ s/(.{45}).*/\1/' ein

# Fix O
#  I=   325 IAn=  8 IBond, RBType=     3  1.00     4  1.00    45  0.10
#  I=   326 IAn=  8 IBond, RBType=     5  1.00     6  1.00    48  0.10    56  0.10
#  I=   328 IAn=  8 IBond, RBType=     9  1.00    10  1.00
#                                                        ^
# sed -r -i '/IAn=  8/ s/(.{57}).*/\1/' ein