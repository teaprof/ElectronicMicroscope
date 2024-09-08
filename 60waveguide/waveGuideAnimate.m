function waveGuideAnimate(sol)
%     if nargin == 0
%         geom = geometry;
%         w = 2*pi*1e+12;
%         m = 0;
%         waveguideSolver = WaveguideSolver(geom, w, m);
%         waveguideSolver = waveguideSolver.solve;
%         sol = waveguideSolver.sol;
%     end
%     save('sol', 'sol');
    load('sol');
    f = figure;
    
    maxr = 2*sol.r(1);
    maxz = 10/sol.kz;
    maxt = 10/sol.w;
    
    rr = linspace(-maxr, maxr, 30);
    zz = linspace(0, maxz, 80);

    subplot(1, 2, 1);
    visPlane1 = mycolorplane([1 0 0], [0 1 0], rr, zz);
    colorbar;
    view(2);
    subplot(1, 2, 2);
    visPlane2 = mycolorplane([1 0 0], [0 1 0], rr, zz);
    colorbar;
    view(2);
    
    tt = linspace(0, maxt, 100);
    rr = [visPlane1.XY(1, :); zeros(1, size(visPlane1.XY, 2)); visPlane1.XY(2, :)];
    for m = 1 : numel(tt)        
        [E, B] = sol.getFieldXYZt_EB(rr, tt(m));
        normE = sqrt(sum(E.^2));
        normB = sqrt(sum(B.^2));
        visPlane1.img.CData(:) = normE(:);
        visPlane2.img.CData(:) = normB(:);
        drawnow;
    end
end