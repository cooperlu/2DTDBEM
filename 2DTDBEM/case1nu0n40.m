clc
clear all
close all
%% INPUT
% Parameters
E = 4/3; % Young's modulus
ro = 2; % density
nu =0; % Possion's Ratio
GE = E/(2*(1+nu)); % Shear Modulus
cs = sqrt(GE/ro);
cp = cs*sqrt((2-2*nu)/(1-2*nu));
at = 0.05; % Time Step Size
nstep = 15/at; % Total Time Step
tend = round(1/at);

% Geometry
% Original Mesh
x = [0:0.05:1,ones([1,19]),1:-0.05:0,zeros([1,19])];
y = [zeros([1,20]),0:0.05:1,ones([1,19]),1:-0.05:0.05];
n = 80;
nn = 2*n;

% Initialize
bcun = zeros([2*n,nstep]); % matrix store solution to Ax = B
bck = zeros([2*n,nstep]); % known boundary condition; NEED MOD
bck(41:2:79,1:tend) = -ones(size(bck(41:2:79,1:tend))); % HARD CODING!!!!!
kode =  zeros([2*n,nstep]); % type of BC
kode(121:160,:) = ones(size(kode(121:160,:)));% HARD CODING!!!!!
udata = bcun; % Matrix Store u data; structure: [sptial x temporal];
pdata = bcun; % Matrix Store p data; structure: [sptial x temporal];
Hdata = zeros([2*n,2*n,nstep]); % Matrix store H^nm; 3rd index = n-m+1; struc: [nn x nn x temporal]
Gdata = zeros([2*n,2*n,nstep]); % Matrix store G^nm; 3rd index = n-m+1; struc: [nn x nn x temporal]

%% CALCULATION
for istep = 1:nstep
    fprintf('Step %5d',istep);
    
    
    if istep == 1
        x(n+1) = x(1);
        y(n+1) = y(1);
        xm = (x(1:end-1)+x(2:end))./2;
        ym = (y(1:end-1)+y(2:end))./2;
        sl = sqrt((xm-x(1:n)).^2 + (ym-y(1:n)).^2);
        etax = (ym - y(1:n))./sl;
        etay = (x(1:n) - xm)./sl;
    end
    %%
    
    for i = 1:n
        for j = 1:n
            if i == j
                [Hws,Gws] = locin11(sl(j),etax(j),etay(j),istep,at,ro,cp,cs);
                if istep  == 1
                    Hws(1,1) = 0.5;
                    Hws(2,2) = 0.5;
                end
                G(2*i-1:2*i,2*j-1:2*j) = Gws(1:2,1:2);
                H(2*i-1:2*i,2*j-1:2*j) = Hws(1:2,1:2);
            else
                dist(i,j) = distc(xm(i),ym(i),x(j),y(j),xm(j),ym(j),sl(j));
                
                [Hws,Gws] = extin11(xm(i),ym(i),x(j),y(j),x(j+1),y(j+1),...
                    sl(j),etax(j),etay(j),at,istep,cp,cs,nu,GE,ro);
                G(2*i-1:2*i,2*j-1:2*j) = Gws(1:2,1:2);
                H(2*i-1:2*i,2*j-1:2*j) = Hws(1:2,1:2);
            end
        end
    end
    
    
    Gdata(1:2*n,1:2*n,istep) = G;
    Hdata(1:2*n,1:2*n,istep) = H;
    if istep == 1
        % A and B is Hnn and Gnn, does not change
        
        A = H;
        B = G;
        for j = 1:nn
            if kode(j,istep) == 1
                A(:,j) = -G(:,j);
                B(:,j) = -H(:,j);
            end
        end
        FP = zeros([2*n,1]);
    else
        % Historical Effect
        FP = zeros([2*n,1]);
        for step = 2:istep
            FP = FP + Gdata(1:2*n,1:2*n,step)*pdata(1:2*n,istep+1-step)-...
                Hdata(1:2*n,1:2*n,step)*udata(1:2*n,istep+1-step);
        end
    end
    
    bcun(1:2*n,istep) = A\(B*bck(1:2*n,istep)+FP);
    % Reorder and save u and p data
    udata(1:2*n,istep) = bcun(1:2*n,istep);
    pdata(1:2*n,istep) = bck(1:2*n,istep);
    for j = 1:nn
        if kode(j,istep) == 1
            udata(j,istep) = bck(j,istep);
            pdata(j,istep) = bcun(j,istep);
        end
    end
end

t = (1:nstep)*at;
plot(t,udata(39,:));
