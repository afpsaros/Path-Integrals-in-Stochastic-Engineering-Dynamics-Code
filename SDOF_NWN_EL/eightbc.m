function res = Xini(ya, yb, m0, c0, k0, e1, e2, p0, q0, r0, S0, tt, BVX)

    res = [ ya(1)
            ya(2)
            ya(3)
            ya(4)
            yb(1) - BVX(1) 
            yb(2) - BVX(2)
            yb(3) - BVX(3)
            yb(4) - BVX(4)];
    
end