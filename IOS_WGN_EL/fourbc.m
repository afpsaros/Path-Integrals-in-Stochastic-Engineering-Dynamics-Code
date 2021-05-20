function res = Xini(ya, yb, m0, c0, k0, e1, e2, S0, tt, BVX)

    res = [ ya(1)
            ya(2)
            yb(1) - BVX(1) 
            yb(2) - BVX(2)];
    
end