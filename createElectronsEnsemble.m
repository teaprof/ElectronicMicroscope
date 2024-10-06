function [r0, v0] = createElectronsEnsemble(v0norm, rmax, zspan, pattern, varargin)
%function [r0, v0] = createElectronsEnsemble(v0norm, rmax, zspan, pattern, [...])
%function [r0, v0] = createElectronsEnsemble(v0norm, rmax, zspan, 'unirandom', N)
%function [r0, v0] = createElectronsEnsemble(v0norm, rmax, zspan, 'normal', N)
%function [r0, v0] = createElectronsEnsemble(v0norm, rmax, zspan, 'regular', N, M)
%function [r0, v0] = createElectronsEnsemble(v0norm, rmax, zspan, 'grid', N, N, n, m)
% Estart - начальная экнергия электронов, эВ
% rmax - радиус облака, м
% zspan - минимальное и макс значение координаты z электронов, м
% pattern - тип распределения, см. исходный код
    z1 = zspan(1);
    z2 = zspan(2);
    switch(pattern)
        case 'unirandom'
            N = varargin{1};
            [r0, v0] = createElectronsRandom(v0norm, rmax, z1, z2, N);
        case 'normal'
            N = varargin{1};
            [r0, v0] = createElectronsNormal(v0norm, rmax, z1, z2, N);
        case 'regular'
            N = varargin{1};
            M = varargin{2};
            [r0, v0] = createElectronsRegular(v0norm, rmax, z1, z2, N, M);
        case 'grid'
            N = varargin{1};
            M = varargin{2};
            n = varargin{3};
            m = varargin{4};
            [r0, v0] = createElectronsGrid(v0norm, rmax, z1, z2, N, M, n, m);
        otherwise
            error("unrecognized pattern")
    end
    idx = findUniqElectrons(r0, 1e-6);
    r0 = r0(:, idx);
    v0 = v0(:, idx);
end

function [r0, v0] = createElectronsRandom(v0norm, r, z1, z2, N)    
    v0 = repmat([0; 0; v0norm], 1, N);
    
    r0 = zeros(3, N);
    for n = 1 : N
        phi = 2*pi*rand;
        rho = r*rand;
        z = z1 + (z2 - z1)*rand;
        r0(:, n) = [rho*cos(phi); rho*sin(phi); z];
    end
end

function [r0, v0] = createElectronsNormal(v0norm, r, z1, z2, N)
    v0 = repmat([0; 0; v0norm], 1, N);
    
    r0 = zeros(3, N);
    meanz = mean([z1 z2]);
    stdz = std2([z1 z2])/2;
    for n = 1 : N
        phi = 2*pi*rand;
        rho = r*rand;
        z = -Inf;
        while(z < z1 || z > z2)
            z = z1 + normrnd(meanz, stdz);
        end
        r0(:, n) = [rho*cos(phi); rho*sin(phi); z];
    end
end


function [r0, v0] = createElectronsRegular(v0norm, r, z1, z2, N, M)
    v0 = repmat([0; 0; v0norm], 1, N*M);
    r0 = zeros(3, N*M);
    for n = 1 : N
        for m = 1 : M
            x = 2*r*((n-0.5)/N - 0.5);
            if(x == 0)
                x = r/N/100; %в нуле EMsolution выдаёт NaN
            end
            y = 0;
            z = z1 + (z2 - z1)*(m-1)/(M-1);
            r0(:, n + (m-1)*N) = [x; y; z];
        end
    end
end

function [r0, v0] = createElectronsGrid(v0norm, r, z1, z2, N, M, n, m)    
    r0 = [];
    for nn = 1 : N
        x = 2*r*((nn-0.5)/N - 0.5);
        if(x == 0)
            x = r/N/100; %в нуле EMsolution выдаёт NaN
        end
        for mm = 1 : m
            y = 0;
            z = z1 + (z2 - z1)*(mm-1)/(m-1);            
            r0(:, end+1) = [x; y; z];
        end
    end
    for mm = 1 : M
        z = z1 + (z2 - z1)*(mm-1)/(M-1);
        for nn = 1 : n
            x = 2*r*((nn-0.5)/n - 0.5);
            if(x == 0)
                x = r/N/100; %в нуле EMsolution выдаёт NaN
            end        
            y = 0;            
            r0(:, end+1) = [x; y; z];
        end
    end
    v0 = repmat([0; 0; v0norm], 1, size(r0, 2));
end

function idx = findUniqElectrons(r, mindr)
    nparticles = size(r, 2);
    idx = zeros(nparticles, 1);
    idxcount = 1;
    idx(1) = 1;
    mindr2 = mindr^2;
    curdist2 = zeros(nparticles, 1);
    
    for n = 2 : nparticles
        curdist2(1:idxcount) = sum((r(:, n) - r(:, idx(1:idxcount))).^2);
        if min(curdist2(1:idxcount)) >= mindr2
            idxcount = idxcount + 1;
            idx(idxcount) = n;
        end
    end
    idx = idx(1:idxcount);
end
