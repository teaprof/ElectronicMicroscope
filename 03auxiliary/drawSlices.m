function drawSlices(f, lb, ub, dim)
    assert(numel(lb) == 3);
    if nargin < 4
        dim = 3;
    end
    count = [100 2 1];
    points = {};
    for k = 1 : dim
        points{k} = linspace(lb(k), ub(k), count(k));
    end
    
    xydim = setdiff([1, 2, 3], dim);
    xdim = xydim(1);
    ydim = xydim(2);
    tdim = dim;
    [xx, yy] = meshgrid(points{xdim}, points{ydim});
    zz = zeros(count(ydim), count(xdim), count(tdim));
    t = points{tdim};
    
    P = ProgressObj(numel(t)*numel(xx));
    for n = 1 : numel(t)
        zzz = xx;
        parfor m = 1 : numel(xx)
            zzz(m) = f([xx(m); yy(m); t(n)]);
            P.increase();
        end
        zz(:, :, n) = zzz;        
    end
    P.done();
        
    zz = zz*1e+3;
    n = 1;
    figure;
    grid on;
    plot(zz(2, end-30:end))
end