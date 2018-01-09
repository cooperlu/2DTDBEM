function [HW,GW] = extin11(xp,yp,x1,y1,x2,y2,sr,eta1,eta2,at,istep,cp,cs,nu,G,ro )

%% Gauss Quadrature Value
Gauss =[1	0.2955242247147529	-0.1488743389816312;
2	0.2955242247147529	0.1488743389816312;
3	0.2692667193099963	-0.4333953941292472;
4	0.2692667193099963	0.4333953941292472;
5	0.2190863625159820	-0.6794095682990244;
6	0.2190863625159820	0.6794095682990244;
7	0.1494513491505806	-0.8650633666889845;
8	0.1494513491505806	0.8650633666889845;
9	0.0666713443086881	-0.9739065285171717;
10	0.0666713443086881	0.9739065285171717;];
gaussw = (Gauss(:,2));
gaussp = (Gauss(:,3));

%%

ax = (x2-x1)/2;
bx = (x2+x1)/2;
ay = (y2-y1)/2;
by = (y2+y1)/2;

GW = zeros(2);
HW = zeros(2);

for i = 1:10
    xco = ax*gaussp(i) + bx - xp;
    yco = ay*gaussp(i) + by - yp;
    ra = sqrt(xco^2 + yco^2);
    if ra > cp*at*istep
    else
        rd1 = xco/ra;
        rd2 = yco/ra;
        rdn = rd1*eta1 + rd2*eta2;
        xja = gaussw(i) * sr;
        [E1,E2,Z1,Z2,Z3] = fundsole(ra,istep,at,cp,cs,nu,G,ro);
        GW(1,1) = GW(1,1) + xja*(E1-E2+2*rd1*rd1*E2);
        GW(1,2) = GW(1,2) + xja*(2*rd1*rd2*E2);
        GW(2,2) = GW(2,2) + xja*(E1-E2+2*rd2*rd2*E2);
        HW(1,1) = HW(1,1) + xja*(rdn*(Z1+2*rd1*rd1*Z2) + eta1*rd1*(Z1+Z3))/ra;
        HW(1,2) = HW(1,2) + xja*(rdn*2*rd1*rd2*Z2 + eta1*rd2*Z1 + eta2*rd1*Z3)/ra;
        HW(2,1) = HW(2,1) + xja*(rdn*2*rd1*rd2*Z2 + eta1*rd2*Z3 + eta2*rd1*Z1)/ra;
        HW(2,2) = HW(2,2) + xja*(rdn*(Z1+2*rd2*rd2*Z2) + eta2*rd2*(Z1+Z3))/ra;    
    end
    
end
GW(2,1) = GW(1,2);
end


