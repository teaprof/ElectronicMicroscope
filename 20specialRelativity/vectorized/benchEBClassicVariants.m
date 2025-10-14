function benchEBClassicVariants
    N = ceil(logspace(1, log10(2000), 10));
    qsource = getElectronCharge;

    legenda = {'Element-wise', 'Vectorized', 'GPU', 'Coder', 'Coder mex'};

    t1 = zeros(numel(N), 1);
    t2 = zeros(numel(N), 1);
    t3 = zeros(numel(N), 1);
    t4 = zeros(numel(N), 1);
    t5 = zeros(numel(N), 1);
    NTimes = 3; %сколько раз вызывать каждый вариант для усреднения времени
    for n = 1 : numel(N)
        r = rand(3, N(n));
        v = rand(3, N(n))*getSpeedOfLight;
    
        for nn = 1 : NTimes
            tic;
            [E1, B1] = getEBClassicElementwise(qsource, r, v);
            t1(n) = t1(n) + toc;

            tic;
            [E2, B2] = getEBClassicVectorized(qsource, r, v);
            t2(n) = t2(n) + toc;    

            tic;
            if canUseGPU
                [E3, B3] = getEBClassicGPU(qsource, r, v);
            else
                E3 = E1;
                B3 = B1;
            end
            t3(n) = t3(n) + toc;    

            tic;
            [E4, B4] = getEBClassicCoder(qsource, r, v);
            t4(n) = t4(n) + toc;    
    
            tic;
            [E5, B5] = getEBClassicCoder_mex(qsource, r, v);
            t5(n) = t5(n) + toc;        
        end
        
        t1 = t1./NTimes;
        t2 = t2./NTimes;
        t3 = t3./NTimes;
        t4 = t4./NTimes;
        t5 = t5./NTimes;

        dE1 = E1 - E2;
        dB1 = B1 - B2;
        dE2 = E1 - E3;
        dB2 = B1 - B3;
        dE3 = E1 - E4;
        dB3 = B1 - B4;
        dE4 = E1 - E5;
        dB4 = B1 - B5;
        assert(max(abs(dE1(:))) < 1e-15);
        assert(max(abs(dB1(:))) < 1e-15);
        assert(max(abs(dE2(:))) < 1e-15);
        assert(max(abs(dB2(:))) < 1e-15);
        assert(max(abs(dE3(:))) < 1e-15);
        assert(max(abs(dB3(:))) < 1e-15);
        assert(max(abs(dE4(:))) < 1e-15);
        assert(max(abs(dB4(:))) < 1e-15);

        fprintf('N particles:  %d\n', N(n));
        fprintf('%10s: %f sec\n', legenda{1}, t1(n));
        fprintf('%10s: %f sec\n', legenda{2}, t2(n));
        if canUseGPU
            fprintf('%10s: %f sec\n', legenda{3}, t3(n));
        end
        fprintf('%10s: %f sec\n', legenda{4}, t4(n));
        fprintf('%10s: %f sec\n', legenda{5}, t5(n));
        
        fprintf('%d out of %d finished\n\n', n, numel(N));        
    end
    fprintf('Tests passed, all Ok\n');

    if ~canUseGPU
        t3 = 0*t3;
        warning('GPU not available');
    end

    f = figure;
    hold on;
    plot(N, t1);
    plot(N, t2);
    plot(N, t3);
    plot(N, t4);
    plot(N, t5);
    legend(legenda);
    grid on;
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    xlabel('n particles');
    ylabel('time, s');
    saveas(f, 'figures/benchEBclassic.fig');
    saveas(f, 'figures/benchEBclassic.png');
end

function [E, B] = getEBClassicElementwise(q, r, v)
    N = size(r, 2);
    E = zeros(3, N);
    B = zeros(3, N);
    for n = 1 : N
        for k = 1 : N
            if k == n
                continue;
            end            
            [Enk, Bnk] = MovingChargeField.getEBclassic(q, r(:, k), v(:, k), r(:, n));
            E(:, n) = E(:, n) + Enk;
            B(:, n) = B(:, n) + Bnk;
        end        
    end
end