function simEMfield
    w = 1e+14;
    geom = generateGeometry('simple');
%     geom = generateGeometry;
    sol = createSolStruct(geom, w);
    
    f = @(x) targetFcn(sol, x);
    [x0, lb, ub] = solToX(sol);
    [x, y] = lsqnonlin(f, x0, lb, ub);
    sol = XtoSol(sol, x);
    
    sol.Ezr{1}(sol.r(1));
    
    [Ez, Er, Ephi, Hz, Hr, Hphi] = getField(1, 0.5*sol.r(1), 0, sol);
    
    f(x)
        
    plotWaveguide(sol);
    
    
    return;
    m = 0;
    a = 1e-4;
    c = getSpeedOfLight;
    assert(w > 2.404*c/a, "Угловая частота должна быть больше критической");
    f = @(x) norm(getEzr(a, w, 1, 0, x(1), m, eps, mju));
    [x, y] = fminsearch(f, 0);
%     nulls = getAllMins(f, [0, 10], 1e+5);
    f(x);
end

function plotWaveguide(sol)
    sol = expandSolution(sol);
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
        Ezr(n) = sol.Ezr{nlayer(n)}(r(n));
        Hzr(n) = sol.Hzr{nlayer(n)}(r(n));
    end    
%     plot(r, Ezr);
    

    phi = linspace(0, 2*pi, 200);
    [x, y] = pol2cart(phi, r);
    for n = 1 : numel(r)
        for m = 1 : numel(phi)
            [Ez(n, m), Er(n, m), Ephi(n, m), Hz(n, m), Hr(n, m), Hphi(n, m)] = getField(nlayer(n), r(n), phi(m), sol);
        end
    end
    plots = {Ez Er Ephi Hz Hr Hphi};
    titles = {'Ez', 'Er', 'Ephi', 'Hz', 'Hr', 'Hphi'};
    for n = 1 : 6
        subplot(2, 3, n);
        mesh(x, y, abs(plots{n}));
        title(titles{n});
    end

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


function sol = createSolStruct(geometry, w)
    c = getSpeedOfLight;
    sol = geometry; %копируем все поля из geometry в sol
    sol.w = w;
    sol.vacuum.k = w/c;
    sol.m = 0;
    for n = 1 : geometry.nlayers
        sol.k2(n) = sol.vacuum.k^2*geometry.eps(n)*geometry.mju(n);
        sol.kxy(n) = 1; %если сюда подставить 0, то аргумент функции Неймана станет 0 и возникнет Inf
        sol.C1(n) = 1; %коэффициент перед функцией Бесселя
        sol.C2(n) = 1; %коэффициент перед функцией Неймана
        sol.C1enabled(n) = 1; %коэффициент перед функцией Бесселя
        sol.C2enabled(n) = 1; %коэффициент перед функцией Неймана
        sol.C1span(n, :) = [-Inf; Inf];
        sol.C2span(n, :) = [-Inf; Inf];
        
        guess = 2.5505 + 1.2474*sol.m;
        besselzero = fzero(@(x) besselj(sol.m, x), guess); %первый нуль функции Бесселя
        sol.kxyspan(n, :) = [0.8*besselzero/sol.r(n); Inf]; 
    end    
    sol.C2enabled(1) = 0; %в центральном слое коэффициент перед функцией Неймана должен быть равен нулю, так как она при r = 0 обращается в бесконечность    
    sol.C1span(1, :) = [1; 1]; %требуем, чтобы этот коэффициент был равен 1, так как он будет задавать амплитуду решения (если учитывать только ГУС, то система уравнений будет однородной => беск много решений)
    sol.C2span(1, :) = [1; 1]; %реально не учитывается, так как C2enabled(1) = 0
end

function y = targetFcn(sol, x)    
    sol = XtoSol(sol, x);
    sol = expandSolution(sol);
    y = 1e+7*getGUSdiscrepancy(sol);
    y = [y; sol.kxy(:)]; %будем также минимизировать kxy - то есть брать решения, соответствующие минимальным корням функции Бесселя
    fprintf('[%f %f %f] : %30.20g\n', x(1), x(2), x(3), sum(y.^2));
end

function [x, lb, ub] = solToX(sol)
    x = zeros(sol.nlayers*3, 1);    
    lb = zeros(sol.nlayers*3, 1);    
    ub = zeros(sol.nlayers*3, 1);    
    for n = 1 : sol.nlayers
        x(n*3-2) = sol.kxy(n);
        x(n*3-1) = sol.C1(n);
        x(n*3  ) = sol.C2(n);
        
        lb(n*3-2) = sol.kxyspan(n, 1);
        ub(n*3-2) = sol.kxyspan(n, 2);
        lb(n*3-1) = sol.C1span(n, 1);
        ub(n*3-1) = sol.C1span(n, 2);
        lb(n*3  ) = sol.C2span(n, 1);
        ub(n*3  ) = sol.C2span(n, 2);
    end
end

function sol = XtoSol(sol, x)
    for n = 1 : sol.nlayers
        sol.kxy(n) = x(n*3-2);
        sol.C1(n) = x(n*3-1);
        sol.C2(n) = x(n*3);
    end
    sol.C1(1) = 1;
    sol = expandSolution(sol);
end

function delta = getGUSdiscrepancy(sol)
    delta = zeros(sol.nlayers*2, 1);
    phi = 0;    
    for n = 1 : sol.nlayers
        r = sol.r(n);
        [Ez1, Er1, Ephi1, Hz1, Hr1, Hphi1] = getField(n, r, phi, sol);
        if n < sol.nlayers
            [Ez2, Er2, Ephi2, Hz2, Hr2, Hphi2] = getField(n+1, r, phi, sol);
        else
            Ez2 = 0;
            Er2 = 0;
            Ephi2 = 0;
            Hz2 = 0;
            Hr2 = 0;
            Hphi2 = 0;
        end
        dE = abs(Ez2 - Ez1) + abs(Ephi2 - Ephi1);
        dH = abs(Hz2 - Hz1) + abs(Hphi2 - Hphi1);
        delta(n*2 - 1) = dE;
        delta(n*2    ) = dH;
    end
end

function res = expandSolution(sol)
    res = sol;
    for n = 1 : sol.nlayers
        res.kz2(n) = sol.vacuum.k^2*sol.eps(n)*sol.mju(n) - sol.kxy(n)^2;
        res.kz(n) = sqrt(res.kz2(n));

        res.Ezr{n} = @(r) getEzr(n, r, sol);
        res.Ezphi{n} = @(phi) getEzphi(n, phi, sol);
        res.Hzr{n} = @(r) 0;
        res.Hzphi{n} = @(phi) 0;
    end
end

function [Ez, Er, Ephi, Hz, Hr, Hphi] = getField(nlayer, r, phi, sol)
    [Er, Ephi, Hr, Hphi] = getTransverse(nlayer, r, phi, sol);
    Ez = sol.Ezr{nlayer}(r)*sol.Ezphi{nlayer}(phi);
    Hz = sol.Hzr{nlayer}(r)*sol.Hzphi{nlayer}(phi);
end

function [Er, Ephi, Hr, Hphi] = getTransverse(nlayer, r, phi, sol)
    Ezr = sol.Ezr{nlayer};
    Ezphi = sol.Ezphi{nlayer};
    Hzr = sol.Hzr{nlayer};
    Hzphi = sol.Hzphi{nlayer};

    kz = sol.kz(nlayer);
    eps = sol.eps(nlayer);
    mju = sol.mju(nlayer);
    w = sol.w;
    coef = 1i*kz/(sol.w^2*eps*mju - kz^2);
    
    %Вычисляем частные производные функций
    %Ez(r, phi) = Ezr(r)*Ezphi(phi)
    %Hz(r, phi) = Hzr(r)*Hzphi(phi)
    dEz_dr = deriv(Ezr, r)*Ezphi(phi);
    dEz_dphi = Ezr(r)*deriv(Ezphi, phi);
    dHz_dr = deriv(Hzr, r)*Hzphi(phi);
    dHz_dphi = Hzr(r)*deriv(Hzphi, phi);
    
    C1 = w*mju/kz;
    %Er = coef*(dEz/dr + w*mju/kz*(1/r*dHz/dphi))
    Er = dEz_dr + C1/r*dHz_dphi;
    Er = coef*Er;
    
    %Ephi = coef*(1/r*dEz/dphi - w*mju/kz*dHz/dr)
    Ephi = 1/r*dEz_dphi - C1*dHz_dr;
    Ephi = coef*Ephi;
    
    C2 = w*eps/kz;
    %Hr = coef*(dHz/dr - w*eps/kz*(1/r*dEz/dphi))    
    Hr = dHz_dr - C2/r*dEz_dphi;
    Hr = coef*Hr;
    
    %Hphi = coef*(1/r*dHz/dphi + w*eps/kz*dEz/dr)
    Hphi = 1/r*dHz_dphi + C2*dEz_dr;
    Hphi = coef*Hphi;
end

function res = deriv(f, x)
    eps = abs(x)*1e-8;
    if(eps == 0)
        eps = 1e-8;
    end
    y1 = f(x + eps);
    y2 = f(x - eps);
    res = (y2 - y1)/2/eps;
end

function Ezphi = getEzphi(nlayer, phi, sol)
    Ezphi = cos(sol.m*phi); %если подставить sin то при sol.m = 0 получаем всегда Ezphi = 0 и Ez = 0
end

function Ezr = getEzr(nlayer, r, sol)
%     phase = w*t + kz*z + m*phi;
%     c = getSpeedOfLight; %m/s
%     k = w/c;
%     kz2 = eps*mju*k^2 - kxy^2;
    m = sol.m;
    Ezr = 0;
    kxy = sol.kxy(nlayer);
    if(sol.C1enabled(nlayer) > 0)
        Ezr = Ezr + sol.C1(nlayer)*besselj(m, kxy*r);
    end
    if(sol.C2enabled(nlayer) > 0)
        Ezr = Ezr + sol.C2(nlayer)*bessely(m, kxy*r);
    end
end
