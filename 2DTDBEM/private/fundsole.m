function [E1,E2,Z1,Z2,Z3] = fundsole(ra,istep,at,cp,cs,nu,G,ro);

E1 = 0;
E2 = 0;
Z1 = 0;
Z2 = 0;
Z3 = 0;
p0 = cp * istep*at/ra;
if p0 < 1
else
    p1 = max([1,cp*(istep-1)*at/ra]);
    p2 = max([1,cp*(istep-2)*at/ra]);
    s0 = max([1,cs*(istep)*at/ra]);
    s1 = max([1,cs*(istep-1)*at/ra]);
    s2 = max([1,cs*(istep-2)*at/ra]);
    cps = 2*(1-nu)/(1-2*nu);
    ps0 = (p0/s0)^2 -1;
    ps1 = (p1/s1)^2 -1;
    ps2 = (p2/s2)^2 -1;
    rais0 = rais(p0) + rais(s0);
    rais1 = rais(p1) + rais(s1);
    rais2 = rais(p2) + rais(s2);
    cte = 1/4/pi/ro/cp/cp;
    E1 = cte * ((acosh(p0) - acosh(p1)) + cps*(acosh(s0)-acosh(s1)));
    E2 = cte * ps0/rais0;
    if p1 > 1
        E2 = E2 - cte*ps1/rais1;
    end
    v1 = istep/3*(rais(p0)*rais0+rais(s0)^2)*ps0/rais0;
    if p1 > 1
        v1 = v1-2*(istep-1)/3*(rais(p1)*rais1+rais(s1)^2)*ps1/rais1;
    end
    if p2 > 1
        v1 = v1+(istep-2)/3*(rais(p2)*rais2+rais(s2)^2)*ps2/rais2;
    end
    v2 = istep*rais(p0)/2 - (istep-1)*rais(p1) + (istep-2)*rais(p2)/2;
    v3 = istep*rais(s0)/2 - (istep-1)*rais(s1) + (istep-2)*rais(s2)/2;
    ctz = 4*G*cte;
    Z1 = ctz*(v1-cps*v3);
    Z2 = ctz*(-2*v1-(v2-cps*v3));
    Z3 = ctz*(v1 - 2*nu/(1-2*nu)*v2);
end
end


function a = rais(x)
a = sqrt(x^2-1)/x;
end