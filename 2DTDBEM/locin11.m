function [HW,GW] = locin11(sr,eta1,eta2,istep,at,ro,cp,cs)
GW = zeros(2);
HW = GW;
time0 = istep*at;
time1 = (istep-1)*at;

p0 = max([1,cp*time0/sr]);
p1 = max([1,cp*time1/sr]);
s0 = max([1,cs*time0/sr]);
s1 = max([1,cs*time1/sr]);
raiz0 = raiz(p0) + p0/s0*raiz(s0);
raiz1 = raiz(p1) + p1/s1*raiz(s1);
act = (time0 * (log(p0+raiz(p0))/p0 + sin(1/p0))-...
    time0 *(log(p1+raiz(p1))/p1 + sin(1/p1)))/cp+...
    (time0 * (log(s0+raiz(s0))/s0 + sin(1/s0))-...
    time0 *(log(s1+raiz(s1))/s1 + sin(1/s1)))/cs;
bct = (time0*sin(1/p0)-time1*sin(1/p1))/cp-...
    (time0*asin(1/s0)-time1*sin(1/s1))/cs;
if p0>1
    bct = bct-time0/cp*(1-(p0/s0)^2)/raiz0;
end
if p1>1
    bct = bct+time1/cp*(1-(p1/s1)^2)/raiz1;
end
glocal11 = (act-bct)/2/pi/ro;
glocal22 = (act+bct)/2/pi/ro;
GW(1,1) = glocal11*eta2^2+glocal22*eta1^2;
GW(1,2) = (glocal22-glocal11)*eta1*eta2;
GW(2,1) = GW(1,2);
GW(2,2) = glocal11*eta1^2+glocal22*eta2^2;


end

function a = raiz(x)
a = sqrt(x^2-1);
end
