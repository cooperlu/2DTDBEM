%%

n = 32;

%%

x(n+1) = x(1);
y(n+1) = y(1);
xm = (x(1:end-1)+x(2:end))./2;
ym = (y(1:end-1)+y(2:end))./2;
sl = sqrt((xm-x(1:n)).^2 + (ym-y(1:n)).^2);
etax = (ym - y(1:n))./sl;
etay = (x(1:n) - xm)./sl;

%%

for i = 1:n
    for j = 1:n
        if i == j;
            locin11(
            if istep  == 1
                HT(1,1) = 0.5;
                HT(2,2) = 0.5;
            end
            G(2*i-1:2*i,2*j-1:2*j) = GT(1:2,1:2);
            H(2*i-1:2*i,2*j-1:2*j) = HT(1:2,1:2);
        else
            dist = distc
            if dist > cp*istip*at
            else
                extin11
                G(2*i-1:2*i,2*j-1:2*j) = GT(1:2,1:2);
                H(2*i-1:2*i,2*j-1:2*j) = HT(1:2,1:2);
            end
        end
        
        
        
    end
end