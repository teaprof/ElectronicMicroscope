function [xx, yy] = optimizator(f, lb, ub, NTrials, Nbest)
%function [xx, yy] = optimizator(f, lb, ub, NTrials, Nbest)
%
% Implement multistart algorithm
    assert(NTrials >= Nbest);
    
    fprintf('Global phase\n');
    seeds = randi(hex2dec('FFFFFFFF'), NTrials, 1);

    P = ProgressObj(NTrials);
    yy = zeros(1, NTrials);
    xx = zeros(numel(lb), NTrials);
    parfor n = 1 : NTrials
        rng(seeds(n))
        x = lb + (ub - lb).*rand(size(lb));
        xx(:, n) = x;
        yy(n) = f(x);
%         opts = optimoptions('fmincon', 'Display', 'none');
%         [xx(:, n) yy(n)] = fmincon(f, x, [], [], [], [], lb, ub, [], opts);
        P.increase(1);
    end
    P.done();
    
    [yy, idx] = sort(yy);
    xx = xx(:, idx);
    
    %Оставляем только лучшие решения
    Nbest = min(NTrials, Nbest);
    yy = yy(1:Nbest);
    xx = xx(:, 1:Nbest);
    return
    
    %"Дожимаем" лучшие решения локальным солвером    
    fprintf('Local phase\n');
    opts = optimoptions('fmincon', 'Algorithm', 'active-set', 'Display', 'final-detailed', 'MaxFunctionEvaluations', 50, 'FunctionTolerance', 1e-5);
    P = ProgressObj(Nbest);
    parfor n = 1 : Nbest
        fprintf('n = %d, seed = %d\n', n, seeds(idx(n)));
        [xxx, yyy] = fmincon(f, xx(:, n), [], [], [], [], lb, ub, [], opts);  
        %К сожалению, вызов fmincon не всегда заканчивается улучшением
        %решения
        if yyy < yy(n)
            xx(:, n) = xxx;
            yy(n) = yyy;
        end
        P.increase(1);
    end
    P.done();
    
    [yy, idx] = sort(yy);    
    xx = xx(:, idx);
end