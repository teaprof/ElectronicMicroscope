%test if cluster works properly.
%see also: makecluster
function clustertest        
    Dim = 10;
    Options = metropolis_get('full');    
    NTasks = numel(Options);
    lb = -10*ones(Dim, 1);
    ub = 10*ones(Dim, 1);
    parfor n = 1 : NTasks
        x0 = lb + (ub-lb).*rand(Dim, 1);            
        val = metropolis(@func, x0, lb, ub, Options{n});
        res(n) = val;
    end      
    for n = 1 : numel(res)
        fprintf('%d: min = %g\n', n, res(n));
    end
end

function res = func(x)
    res = norm(x) + 3*sin(abs(norm(x)));
end