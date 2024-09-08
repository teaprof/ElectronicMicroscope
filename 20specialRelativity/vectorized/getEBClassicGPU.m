function [E, B] = getEBClassicGPU(q, r, v)
    if nargin == 0
        N = 10;
        r = rand(3, N);
        v = rand(3, N);
        q = getElectronCharge;
    end
    [E, B] = getEBClassicVectorized(q, gpuArray(r), gpuArray(v));
    E = gather(E);
    B = gather(B);
end

