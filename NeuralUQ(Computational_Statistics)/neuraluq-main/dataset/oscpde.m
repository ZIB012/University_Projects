function [c,f,s] = oscpde(x,t,u,dudx)
    nu = 0.5;
    
    c = 1;
    f = nu*dudx;
    s = -u*dudx;
end
