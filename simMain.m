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
%     electronsPerParticle = 1e+4;    
%     save('simMain');
%     return;
    load('simMain');
    vphase = (0.95 + 0.05*rand)*getSpeedOfLight;
    sol.kz = sol.w/vphase;
    lambda = 2*pi/sol.kz;
    
    EMScale = 1e+7; %по умолчанию амплитуда Ez равна 1 В/м, этот множитель масштабирует решение
    waveguideSolver.sol.C1 = waveguideSolver.sol.C1*EMScale;
    waveguideSolver.sol.C2 = waveguideSolver.sol.C2*EMScale;
    
    %[r0, v0] = createElectronsEnsemble(1, geom.r(1)*2, [0 lambda], 'unirandom', 1e+2);
    [r0, v0] = createElectronsEnsemble(1, geom.r(1)*2, [0 lambda], 'regular', 19, 299);
    [r0, v0] = createElectronsEnsemble(1, geom.r(1)*2, [0 lambda], 'regular', 19, 19);
    %[r0, v0] = createElectronsEnsemble(1, geom.r(1)*9.3, [0 lambda/4], 'grid', 5, 5, 19, 19);
    %[r0, v0] = createElectronsEnsemble(1, geom.r(1), [0 1e-4], 'normal', 1e+2);
    
    tspan = [0 100]*1e-12;   
    
    tic
    v0 = v0*vphase;
    traj = simElectrons(r0, v0, waveguideSolver.sol, tspan, electronsPerParticle);
    toc
    
    
    save('simMainTraj','traj');
    animateFlight2(waveguideSolver.sol, traj);
end
