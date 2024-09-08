function sol = simMain
%     geom = geometry('simple');
%     geom = geometry;
%     w = 2*pi*1e+12;
%     m = 0;
%     waveguideSolver = WaveguideSolver(geom, w, m);
%     waveguideSolver = waveguideSolver.solve;
%     for n = 1 : min(numel(waveguideSolver.solutions), 10)
%         waveguidePlot(waveguideSolver.solutions{n});
%     end
%     sol = waveguideSolver.sol;
%     save('simMain');
%     return;
    load('simMain');
    vphase = (0.8 + 0.1*rand)*getSpeedOfLight;
    sol.kz =  sol.w/vphase;
    lambda = 2*pi/sol.kz;
    
    EMScale = 1e+7; %по умолчанию амплитуда Ez равна 1 В/м, этот множитель масштабирует решение
    waveguideSolver.sol.C1 = waveguideSolver.sol.C1*EMScale;
    waveguideSolver.sol.C2 = waveguideSolver.sol.C2*EMScale;
    
    [r0, v0] = createElectronsEnsemble(geom.r(1)*2, [0 lambda], 'random', 1e+3);
    [r0, v0] = createElectronsEnsemble(geom.r(1)*2, [0 lambda], 'regular', 19, 299);
    [r0, v0] = createElectronsEnsemble(geom.r(1)*9.3, [0 lambda/4], 'grid', 5, 5, 19, 19);
    
    tspan = [0 100]*1e-12;
    
    tic
    traj = simElectorns(r0, v0, waveguideSolver.sol, tspan);
    toc
    
    save('simMainTraj','traj');
    animateFlight2(waveguideSolver.sol, traj);
end


