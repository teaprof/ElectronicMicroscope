function drawEM
%function drawEM
%   Visualizes first n modes of the waveguide solutions

    %Initialize geometry
    %geom = geometry('simple');
    geom = geometry('waveguide');
    
    % Show some EM solutions
    w = 2*pi*1e+12;
    m = 0;
    waveguideSolver = WaveguideSolver(geom, w, m);
    waveguideSolver = waveguideSolver.solve(3);
    for n = 1 : numel(waveguideSolver.solutions)
        filename_base = sprintf('figures/TM_%02d_%02d_mode', waveguideSolver.solutions{n}.m, n);
        waveguidePlot(waveguideSolver.solutions{n}, filename_base);
    end
end
