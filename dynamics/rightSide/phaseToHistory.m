function [xhist, vhist] = phaseToHistory(xhistory, particlen)
    N = size(xhistory, 1)/6;
    idx = (particlen-1)*3 + 1;
    xhist = xhistory(idx:idx+2, :);
    idx = 3*N + idx;
    vhist = xhistory(idx:idx+2, :);
end
