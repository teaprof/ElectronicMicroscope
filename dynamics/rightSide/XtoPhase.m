function [r, v] = XtoPhase(x, N)
    if nargin == 1
        N = numel(x)/6;
        assert(N*6 == numel(x));
    end
    r = reshape(x(1:3*N), 3, N);
    v = reshape(x(3*N+1:end), 3, N);
end
