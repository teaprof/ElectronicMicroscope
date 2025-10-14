function sol = simMain
    %Initialize geometry
    %geom = geometry('simple');
    geom = geometry('waveguide');
    
    % Show some EM solutions
    w = 2*pi*1e+12;
    m = 0;
    waveguideSolver = WaveguideSolver(geom, w, m);
    waveguideSolver = waveguideSolver.solve;
    sol = waveguideSolver.sol;

    %Simulate flight
    electronsPerParticle = 1e+4;
    %vphase = (0.95 + 0.05*rand)*getSpeedOfLight;
    %sol.kz = sol.w/vphase;
    vphase = sol.w/sol.kz;
    vgroup = getSpeedOfLight*(getSpeedOfLight/vphase);
    lambda = 2*pi/sol.kz;
    
    EMScale = 1e+10; %по умолчанию амплитуда Ez равна 1 В/м, этот множитель масштабирует решение
    sol.C1 = sol.C1*EMScale;
    sol.C2 = sol.C2*EMScale;
    
    %[r0, v0] = createElectronsEnsemble(1, geom.r(1)*2, [0 lambda], 'unirandom', 1e+2);
    %[r0, v0] = createElectronsEnsemble(1, geom.r(1)*0.5, [0 lambda], 'regular', 19, 299);
    %[r0, v0] = createElectronsEnsemble(1, geom.r(1)*0.5, [0 lambda], 'regular', 19, 19);
    %[r0, v0] = createElectronsEnsemble(1, geom.r(1)*9.3, [0 lambda/4], 'grid', 5, 5, 19, 19);
    [r0, v0] = createElectronsEnsemble(1, geom.r(1)*0.5, [0 0.2*lambda], 'normal', 4e+1);
    v0 = v0*vgroup*0.95;
    
    tspan = [0 100]*1e-12;
    
    tic    
    traj = simElectrons(r0, v0, sol, tspan, electronsPerParticle);
    toc    
    
    animateFlight2(sol, traj, 'rmax', max(sol.r)*1.1, 'zhalfwidth', lambda, 'gridSize', [40 30], 'fileName', '90figures/simMain.avi');
end
