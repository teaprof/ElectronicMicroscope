function [X, Y, Z, VX, VY, VZ] = solToPhase(odesol)
%See also: XtoPhase, PhaseToX
    x = odesol.y;
    N = size(x, 1)/6;
    idx = 1:3:3*N;
    X = x(idx, :);
    Y = x(idx+1, :);
    Z = x(idx+2, :);
    idx = 3*N + idx;
    VX = x(idx, :);
    VY = x(idx+1, :);
    VZ = x(idx+2, :);
end

