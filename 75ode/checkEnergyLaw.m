function checkEnergyLaw(traj)
    %%Рисуем энергию от времени
    N = traj.Nparticles;
    for n = 1 : numel(traj.t)
        curr = reshape(traj.rv(1:3*N, n), 3, N);
        curv = reshape(traj.rv(3*N+1:6*N, n), 3, N);
        E(n) = sum(getEnergy(curr, curv));
    end    
    figure;
    plot(traj.t, E, '.-');
    fprintf('dE = %g%%\n', (max(E) - min(E))/mean(E)*100);              
end

function E = getEnergy(r, v)
    qsource = getElectronCharge;
    m = getElectronMass;
    N = size(r, 2);
    E = zeros(N, 1);
    for n = 1 : N
        phi = 0;
        for k = 1 : n - 1 %каждая пара частиц (n, k) должна учитываться только один раз
            phi = phi + MovingChargeField.getPotentialClassic(qsource, r(:, k), r(:, n));
        end
        T = m*(norm(v(:, n)))^2/2;
        E(n) = qsource*phi + T;
    end
end
