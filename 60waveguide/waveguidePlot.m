function WaveguidePlot(sol) %рисовалка
    if nargin == 0
        geom = geometry;
        w = 2*pi*1e+12;
        m = 0;
        waveguideSolver = WaveguideSolver(geom, w, m);
        waveguideSolver = waveguideSolver.solve;
        sol = waveguideSolver.sol;
    end
%     WaveguidePlot3D(sol);
%     WaveguidePlot2D(sol);
    WaveguidePlot2D_2(sol);
end

function WaveguidePlot2D(sol) %рисовалка
    %     figure;
    npointsPerLayer = 300;
    r = zeros(sol.nlayers*npointsPerLayer, 1);
    nlayer = zeros(sol.nlayers*npointsPerLayer, 1);
    rmin = min(sol.r(1))*1e-3;
    for n = 1 : sol.nlayers
        rmax = sol.r(n);
        rr = linspace(rmin, rmax, npointsPerLayer);
        r((n-1)*npointsPerLayer + 1: n*npointsPerLayer) = rr;
        nlayer((n-1)*npointsPerLayer + 1: n*npointsPerLayer) = n;
        rmin = rmax;
    end

    for n = 1 : numel(r)        
        Ezr(n) = sol.getEzr(nlayer(n), r(n));
        [Ez(n), Er(n), Ephi(n), Hz(n), Hr(n), Hphi(n)] = sol.getComplexAmplitude_EH(nlayer(n), r(n), 0);
    end
    figure;
    subplot(2, 1, 1);
    plot(r, Ezr);
    xlabel('r, m');
    ylabel('V/m');
    title('E_z(r)');
    grid('on');
    subplot(2, 1, 2);
    plot(r, imag(Hphi));
    title('H_{\phi}(r)');
    xlabel('r, m');
    ylabel('A/m');
    grid('on');
    
    
    
    subplot(2, 1, 1);
    hold('on');
    for n = 1 : sol.nlayers
        plot([sol.r(n), sol.r(n)], [min(Ezr), max(Ezr)]);
    end
end

function WaveguidePlot2D_2(sol) %рисовалка
    %     figure;
    npointsPerLayer = 300;
    r = zeros(sol.nlayers*npointsPerLayer, 1);
    nlayer = zeros(sol.nlayers*npointsPerLayer, 1);
    rmin = min(sol.r(1))*1e-3;
    for n = 1 : sol.nlayers
        rmax = sol.r(n);
        rr = linspace(rmin, rmax, npointsPerLayer);
        r((n-1)*npointsPerLayer + 1: n*npointsPerLayer) = rr;
        nlayer((n-1)*npointsPerLayer + 1: n*npointsPerLayer) = n;
        rmin = rmax;
    end

    for n = 1 : numel(r)        
        [Ez(n), Er(n), Ephi(n), Hz(n), Hr(n), Hphi(n)] = sol.getComplexAmplitude_EH(nlayer(n), r(n), 0);
    end
    figure;
    subplot(2, 1, 1);
    plot(r, Ez);
    xlabel('r, m');
    ylabel('V/m');
    title('E_z(r)');
    grid('on');
    subplot(2, 1, 2);
    plot(r, imag(Er));
    title('E_{r}(r)');
    xlabel('r, m');
    ylabel('V/m');
    grid('on');
end


function WaveguidePlot3D(sol) %рисовалка
    %     figure;
    npointsPerLayer = 20;
    r = zeros(sol.nlayers*npointsPerLayer, 1);
    nlayer = zeros(sol.nlayers*npointsPerLayer, 1);
    rmin = 0;
    for n = 1 : sol.nlayers
        rmax = sol.r(n);
        rr = linspace(rmin, rmax, npointsPerLayer);
        r((n-1)*npointsPerLayer + 1: n*npointsPerLayer) = rr;
        nlayer((n-1)*npointsPerLayer + 1: n*npointsPerLayer) = n;
        rmin = rmax;
    end

    for n = 1 : numel(r)
        Ezr(n) = sol.getEzr(nlayer(n), r(n));
        Hzr(n) = sol.getHzr(nlayer(n), r(n));
    end
    %     plot(r, Ezr);


    phi = linspace(0, 2*pi, 200);
    [x, y] = pol2cart(phi, r);
    for n = 1 : numel(r)
        for m = 1 : numel(phi)
            [Ez(n, m), Er(n, m), Ephi(n, m), Hz(n, m), Hr(n, m), Hphi(n, m)] = sol.getComplexAmplitude_EH(nlayer(n), r(n), phi(m));
        end
    end
    plots = {Ez Er Ephi Hz Hr Hphi};
    figure;
    titles = {'Ez', 'Er', 'Ephi', 'Hz', 'Hr', 'Hphi'};
    for n = 1 : 6
        subplot(2, 3, n);
        mesh(x, y, abs(plots{n}));
        title(titles{n});
    end
    return

    figure;
    Ex = abs(Er).*cos(phi);
    Ey = abs(Er).*sin(phi);
    quiver(x, y, Ex, Ey);
    axis equal;

    figure;
    Hx = -abs(Hphi).*sin(phi);
    Hy = abs(Hphi).*cos(phi);
    quiver(x, y, Hx, Hy);
    axis equal;
end
