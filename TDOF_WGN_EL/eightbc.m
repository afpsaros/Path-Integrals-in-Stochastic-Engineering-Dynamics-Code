function res = Xini(ya, yb, m, c, k, e1, e2, S0, tt, BVX1, BVX2)

    res = [ ya(1)
            ya(2)
            yb(1) - BVX1(1) 
            yb(2) - BVX2(1)
            ya(5)
            ya(6)
            yb(5) - BVX1(2) 
            yb(6) - BVX2(2)];
    
end