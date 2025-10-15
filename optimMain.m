function optimMain
    rng(1);
    Clock = clock;
    pathToSave = sprintf('optimresults/run%4d-%02d-%02d-%02d-%02d-%02d', floor(Clock));
    mkdir(pathToSave);
    
%     assert(vstart < getSpeedOfLight/1.05);

    range = ones(getDim, 2)*NaN; %пока инициализируем NaN'ами
    range(getDim('r1'), :) = [1e-3 2e-3]; %м
    range(getDim('r2_minus_r1'), :) = [1e-4 1e-3]; %м
    range(getDim('EMamplitude'), :) = [1e+5 1e+9]; %В/м
    range(getDim('EMphase'), :) = pi + [pi/2-0.1 pi/2+-0.1]; %рад
    range(getDim('EMphase'), :) = [0, 2*pi]; %рад
    range(getDim('EStart'), :) = [0.6 200]*1e+6; %эВ
    range(getDim('tmax'), :) = [30e-12 30e-12]; %с
    assert(prod(prod(~isnan(range))) == 1); %если ничего не забыли проинициализировать, то NaN быть не должно
    lb = range(:, 1).*getScale;
    ub = range(:, 2).*getScale;
    
    electronsPerParticle = 1e+2; % этот параметр не нужен, так как соотношение заряда и массы не зависит от числа электронов, формирующих макрочастицу
    
    %Создаём ансамбль частиц, распределённых в цилиндре радиусом 1 и
    %высотой 1. Этот ансамбль будет потом в соответствующих пропорциях
    %сжиматься, чтобы поместиться в волновод. То же самое касается
    %начальной скорости.
    %Если генерить каждый раз новый ансамбль, то в расчётах значения целевой функции будет
    %статистическая погрешность, которая не даст нормально минимизировать
    %её значение традиционными методами.
    % [r0, v0] = createElectronsEnsemble(1, 0.9, [0 1], 'grid', 5, 5, 19, 19);
    %[r0, v0] = createElectronsEnsemble(1, 0.9, [0 1e-4], 'regular', 5, 5);
    [r0, v0] = createElectronsEnsemble(1, 0.5, [0 0.4], 'unirandom', 2e+2);
    %[r0, v0] = createElectronsEnsemble(1, 0.5, [-0.5 0.5], 'normal', 1e+2); %todo: magic values 1->waveguide inner radius, 1e-4 -> waveguide len
    
    NParticles = size(r0, 2);
    
    charge = NParticles*electronsPerParticle*getElectronCharge;
    fprintf('Charge %f pC\n', charge*1e+12);
    
   
    %Глобальный поиск экстремума
    NTrials = 100; %количество статистических испытаний для глобального поиска
    Nbest = min(NTrials, 8); %количество решений, которые будут дожиматься локальным оптимизатором
    f = @(x)targetFcnExtended(x, r0, v0, electronsPerParticle, lb, ub);
    [xx, yy] = optimizator(f, lb, ub, NTrials, Nbest);
    %drawSlices(f, lb, ub, {'vphase', 'E', 'phase0', 'tspan'}) %тут можно вызвать другой визуализатор todo: переписать функцию для многомерного использования, а не только для dim=3    
    
    %Печатаем таблицу лучших решений 
    for n = 1 : Nbest
        fprintf('#%d: x = [',  n);
        for k = 1 : size(xx, 1)
            fprintf(' %f', xx(k, n));
        end
        fprintf('], y = %f mm\n', yy(n)*1e+3);
    end    
               
    P = ProgressObj(Nbest);
    for n = 1:Nbest %sometimes parfor here doesn't work
        suffix = sprintf('%02d', n);
        OutputSolution(xx(:, n), r0, v0, electronsPerParticle, pathToSave, suffix);
        P.increase(1);
    end
    P.done();
end

function OutputSolution(x, r0, v0, electronsPerParticle, pathToSave, suffix)
    %Пересчитываем траекторию с большей детализацией
    fprintf('Reconstructing optimal config track\n');
    [sol, r, diam, tmax, ~] = XtoStruct(x, r0, v0);
    traj = simflight(sol, r, diam, tmax, electronsPerParticle);    


    y0 = getEnsembleLength(traj.Z(:, 1));
    y = targetFcn(traj);
    fprintf('Initial length %f mm\n', y0*1e+3);
    fprintf('Final length %f mm\n', y*1e+3);
    fprintf('Compression level %f\n', y/y0);
    
    diam = zeros(numel(traj.t), 1);
    for n = 1 : numel(traj.t)
        diam(n) = getEnsembleLength(traj.Z(:, n));
    end
    f = figure;
    plot(traj.t*1e+12, diam*1e+3);
    xlabel('t, ps');
    ylabel('length, mm');
    grid on;
    saveas(f, [pathToSave, '/ensembleLength' suffix]);
    
    save([pathToSave, '/optimMainRes' suffix]); %, 'x', 'y', 'sol', 'traj');
    animateFlight(sol, traj, ...
        'rmax', max(sol.r)*1.1, ...
        'zhalfwidth', max(diam), ...
        'gridSize', [40 30], ...
        'figSize', [800 600], ...
        'aviFileName', [pathToSave '/Fig' suffix '.avi'], ...
        'gifFileName', [pathToSave '/Fig' suffix '.gif']);
    
    %Сохраняем результаты в виде текста
    Estart = getValue(x, 'EStart')/getValue(getScale(), 'EStart');
    f = fopen([pathToSave, '/params', suffix, '.txt'], 'w+');    
    fprintf(f, 'E0 = %g MeV\n', Estart*1e-6);
    fprintf(f, 'V0 = %g *c\n', energyToSpeed(Estart)/getSpeedOfLight);
    fprintf(f, 'kz = %g 1/m\n', sol.kz);
    fprintf(f, 'Ez = %g V/m\n', sol.C1(1));
    fprintf(f, 'phase = %f rad\n', sol.phase);
    fprintf(f, 'Initial length %f mm\n', y0*1e+3);
    fprintf(f, 'Final length %f mm\n', y*1e+3);
    fprintf(f, 'Compression level %f\n', y/y0);
    fclose(f);    
end

function idx = getDim(name)
    if nargin == 0
        %Общее количество элементов в фазовом пространстве
        idx = 6;
    else
        %Номер размерности для заданного имени
        switch(name)
            case 'r1'
                idx = 1;
            case 'r2_minus_r1'
                idx = 2;
            case 'EMamplitude'
                idx = 3;
            case 'EMphase'
                idx = 4;
            case 'EStart'
                idx = 5;
            case 'tmax'
                idx = 6;
            otherwise
                error('unknown name %s', name);
        end
    end
end

function val = getValue(x, name)
    val = x(getDim(name));
end
    
function scale = getScale
%используется, чтобы привести все переменные, по которым проводится
%оптимизация, к одному порядку величины
    scale = zeros(getDim, 1);
    scale(getDim('r1')) = 1e+3;
    scale(getDim('r2_minus_r1')) = 1e+3;
    scale(getDim('EMamplitude')) = 1e-4;
    scale(getDim('EMphase')) = 1;
    scale(getDim('EStart')) = 1e-6;
    scale(getDim('tmax')) = 1e+12;
    assert(prod(scale) ~= 0); %если ничего не забыли проинициализировать, то нулей быть не должно
end

function [EMsolution, r, v, tmax, x] = XtoStruct(x, r0, v0)
    x = x./getScale;
    r1 =  getValue(x, 'r1'); %внутренний диаметр волновода
    r2_minus_r1 =  getValue(x, 'r2_minus_r1'); %толщина диэлектрической стенки волновода
    EMamplitude = getValue(x, 'EMamplitude'); %по умолчанию амплитуда Ez равна 1 В/м, этот множитель масштабирует решение в волноводе
    EMphase = getValue(x, 'EMphase'); %фаза волны
    Estart = getValue(x, 'EStart'); %Стартовая энергия электронов
    tmax = getValue(x, 'tmax'); %длительность полёта
    
    r2 = r1 + r2_minus_r1;
    v0norm = energyToSpeed(Estart);
    
    
    %Находим нужную моду решения в волноводе
    geom = geometry('none');
    geom = geom.addLayer(1, 1, r1); %todo: make as parameter
    geom = geom.addLayer(5.5, 1, r2); %todo: make as parameter
    w = 2*pi*1e+12; % Hz  %todo: make as parameter
    m = 0;
    waveguideSolver = WaveguideSolver(geom, w, m);
    
    waveguideSolver = waveguideSolver.solve;
%    for n = 1 : min(numel(waveguideSolver.solutions), 10)
%         waveguidePlot(waveguideSolver.solutions{n});
%     end
    EMsolution = waveguideSolver.sol;
    
    
    EMsolution.C1 = EMsolution.C1*EMamplitude;
    EMsolution.C2 = EMsolution.C2*EMamplitude;
    EMsolution.phase = EMphase;
    
    r = r0;
    r(1, :) = r(1, :)*geom.r(1);
    r(2, :) = r(2, :)*geom.r(1);
    r(3, :) = r(3, :)*1e-4; %todo: make as parameter

    v = v0*v0norm;
end

function x = structToX(r1, dr, EMamplitude, EMphase, tmax)
    x = zeros(3, 1);
    x(getDim('r1')) = r1;
    x(getDim('r2_minus_r1')) = dr;
    x(getDim('EMamplitude')) = EMamplitude;
    x(getDim('EMphase')) = EMphase;
    x(getDim('tmax')) = tmax;
    x = x.*getScale;
end

%function [v, EMsolution, traj] = targetFcnExtended(x, r0, v0, electronsPerParticle)
function res = targetFcnExtended(x, r0, v0, electronsPerParticle, lb, ub)
%fmincon function can violate constrains at intermediate iterations. We
%expand our target function beyond [lb, ub] 
    [EMsolution, r0, v0, tmax] = XtoStruct(x, r0, v0);
    %Check violation of boundaries for tmax only. If tmax < 0 then ode45
    %solver produces unfeasible points
    x(x < lb) = lb(x < lb);
    x(x > ub) = ub(x > ub);
    if tmax == 0
        res = getEnsembleLength(r0);
        return;
    end
    tspan = [0 tmax];
    traj = simElectrons(r0, v0, EMsolution, tspan, electronsPerParticle);
    res = targetFcn(traj);
end
   
   
%function [v, EMsolution, traj] = targetFcn(EMsolution, r0, v0, tmax, electronsPerParticle)    
function res = targetFcn(traj)
    % вычисляем длину (диаметр) ансамбля как функцию времени
    diam = zeros(numel(traj.t), 1);
    for n = 1 : numel(traj.t)
        diam(n) = getEnsembleLength(traj.Z(:, n));
    end
    %Результат - минимальная длина ансамбля
    res = min(diam);
    
    %Старый способ:
    %Результат - это средня длина ансамбля за последнюю p-часть  времени
    %v = griddedInterpolant(traj.t, diam);
    %p = 0.2; % какую часть времени, начиная с хвоста, учитывать при вычислении среднего
    %t = linspace(tspan(1) + diff(tspan)*(1-p), tspan(2), 1000);
    %res = max(len(t));
end

function traj = simflight(EMsolution, r0, v0, tmax, electronsPerParticle)    
    tspan = [0 tmax];
    traj = simElectrons(r0, v0, EMsolution, tspan, electronsPerParticle);
end

function v = getEnsembleLength(Z)    
    % fast version:
    v = max(Z) - min(Z);
    % robust version:
    %v = 2.0*std2(Z);
end
