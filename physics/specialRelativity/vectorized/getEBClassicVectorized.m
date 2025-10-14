function [E, B] = getEBClassicVectorized(q, r, v)
    if nargin == 0
        N = 10;
        r = rand(3, N);
        v = rand(3, N);
        q = getElectronCharge;
    end
    eps0 = getEps0;
    mju0 = getMju0;
    
    N = size(r, 2);
    
    dr(1, :, :) = r(1, :) - r(1, :)';
    dr(2, :, :) = r(2, :) - r(2, :)';
    dr(3, :, :) = r(3, :) - r(3, :)';
    
    norm_dr3 = squeeze((sum(dr.^2, 1)).^1.5);
    
    
    Etemp = zeros(3, N, N);
    Etemp(1, :, :) = squeeze(dr(1, :, :))./norm_dr3;
    Etemp(2, :, :) = squeeze(dr(2, :, :))./norm_dr3;
    Etemp(3, :, :) = squeeze(dr(3, :, :))./norm_dr3;
    for n = 1 : N
        Etemp(1, n, n) = 0;
        Etemp(2, n, n) = 0;
        Etemp(3, n, n) = 0;
    end
    E = 1/4/pi/eps0*squeeze(sum(Etemp, 2))*q;
    
    vv = repmat(v, 1, 1, N);
    
    Cross = cross(vv, dr);    
    Btemp = zeros(3, N, N);
    Btemp(1, :, :) = squeeze(Cross(1, :, :))./norm_dr3;
    Btemp(2, :, :) = squeeze(Cross(2, :, :))./norm_dr3;
    Btemp(3, :, :) = squeeze(Cross(3, :, :))./norm_dr3;
    for n = 1 : N
        Btemp(1, n, n) = 0;
        Btemp(2, n, n) = 0;
        Btemp(3, n, n) = 0;
    end
    B = mju0/4/pi*squeeze(sum(Btemp, 2))*q;
end

