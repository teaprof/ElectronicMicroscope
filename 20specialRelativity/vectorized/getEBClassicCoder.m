function [E, B] = getEBClassicCoder(q, r, v)
    if nargin == 0
        N = 1e+3;
        r = rand(3, N);
        v = rand(3, N);
        q = getElectronCharge;
    end
    eps0 = getEps0;
    mju0 = getMju0;
    
    N = size(r, 2);
    E = zeros(3, N);
    B = zeros(3, N);
    Emult = 1/4/pi/eps0*q;
    Bmult = mju0/4/pi*q;
    for n = 1 : N
        for k = n + 1 : N
            dr = r(:, n) - r(:, k);
%             norm_dr3 = sum(dr.^2)^1.5;
            norm_dr3 = norm(dr)^3; %fastest
%             norm_dr3 = sqrt(sum(dr.^2))^3;
            EE = dr/norm_dr3;
            E(:, n) = E(:, n) + Emult*EE;
            E(:, k) = E(:, k) - Emult*EE;
            
            c = cross(v(:, k), dr);
            B(:, n) = B(:, n) + Bmult*c/norm_dr3;
            
            c = -cross(v(:, n), dr);
            B(:, k) = B(:, k) + Bmult*c/norm_dr3;
            
        end        
    end
end

