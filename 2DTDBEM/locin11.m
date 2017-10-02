function [GW,HW] = locin11(sr,eta1,eta2,istep,at,ro,cp,cs)
HW = zeros(2);
time0 = istep*at;
time1 = (istep-1)*at;

p0 = max([1,cp*time0/sr]);
p1 = max([1,cp*time1/sr]);
s0 = max([1,cs*time0/sr]);
s1 = max([1,cs*time1/sr]);
raiz0 = raiz(p0) + p0/s0*raiz(s0);
raiz1 = raiz(p1) + p1/s1*raiz(s1);


end

function a = raiz(x)
a = sqrt(x^2-1);
end
