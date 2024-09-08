function animateFlight2(EMfieldSolution, traj, varargin)
    %Парсим параметры
%     load;
    parser = inputParser;
    parser.addParameter('zhalfwidth', 1e-5, @isnumeric);
    parser.addParameter('rmax', 3e-3, @isnumeric);
    parser.addParameter('gridSize', [40 30], @(x) isnumeric(x) && (numel(x) == 2));
    parser.addParameter('figSize', [800 600], @(x) isnumeric(x) && (numel(x) == 2));
    parser.addParameter('fileName', 'WangFig1.avi', @isstr);
    parser.parse(varargin{:});
    params = parser.Results;
    
   
    gridSize = params.gridSize; %[Y, X]
    gridSizeZ = gridSize(1);
    gridSizeX = gridSize(2);
    gridSizeY = gridSize(2);
    
%     params.rmax = 0.01;
%     params.zhalfwidth = 0.01;
    
    %Создаём новое окно
    fig = figure;
    fig.Position = [100 100 params.figSize(1) params.figSize(2)];
    ax1 = subplot(2, 2, 2);
    ax2 = subplot(2, 2, 4);
    axside1 = subplot(2, 2, 1);
    axside2 = subplot(2, 2, 3);
    hold(ax1, 'on');
    hold(ax2, 'on');    
    Column1Width = 0.3;
    Column2Width = 0.6;
    assert(Column1Width + Column2Width < 0.9);
    ax1.OuterPosition(1) = Column1Width;
    ax2.OuterPosition(1) = Column1Width;
    ax1.OuterPosition(3) = Column2Width;
    ax2.OuterPosition(3) = Column2Width;
    axside1.OuterPosition(1) = 0.01;
    axside1.OuterPosition(3) = Column1Width - axside1.Position(1);
    axside2.OuterPosition(1) = 0.01;
    axside2.OuterPosition(3) = Column1Width - axside2.Position(1);
    hold(axside1, 'on');
    hold(axside2, 'on');
    xlabel(ax1, 'z, m');
    ylabel(ax1, 'x, m');
%     xlabel(ax2, 'z, m');
%     ylabel(ax2, 'x, m');
    xlabel(axside1, 'y, m');
    ylabel(axside1, 'x, m');
    xlabel(axside2, 'x, m');
    ylabel(axside2, 'y, m');
    zlabel(axside2, 'z, m');
    view(axside2, 3);
    
    colormap(ax1, MyColorMap);
    colormap(ax2, WongColorMap);
    colormap(axside1, MyColorMap);
    colormap(axside2, WongColorMap);

    %Достаём нужные переменные
    [X, Y, Z] = deal(traj.X, traj.Y, traj.Z);
    Nparticles = traj.Nparticles;
    t = traj.t;
    
    %центр ансамбля - используется для наведения камеры
    CX = griddedInterpolant(t, mean(X, 1));
    CY = griddedInterpolant(t, mean(Y, 1));
    CZ = griddedInterpolant(t, mean(Z, 1));
    params.visualXRange = [-params.rmax; params.rmax]; %[m] поле зрения вокруг центральной точки по X (вертикальная ось)
    params.visualYRange = [-params.rmax; params.rmax]; %[m] поле зрения вокруг центральной точки по X (вертикальная ось)
    params.visualZRange = [-params.zhalfwidth; params.zhalfwidth]; %[m] поле зрения вокруг центральной точки по Z (горизонтальная ось)
    
    
    %Плоскости для визуализации поля
    gridX = linspace(params.visualXRange(1), params.visualXRange(2), params.gridSize(1));
    gridY = linspace(params.visualYRange(1), params.visualYRange(2), params.gridSize(1));
    gridZ = linspace(params.visualZRange(1), params.visualZRange(2), params.gridSize(2));
    [meshXYx, meshXYy] = meshgrid(gridX, gridY);
    [meshYZy, meshYZz] = meshgrid(gridY, gridZ);
    [meshZXx, meshZXz] = meshgrid(gridX, gridZ);
    
    pXY = PhysicalPlane(meshXYx, meshXYy, 0*meshXYx);
    pYZ = PhysicalPlane(0*meshYZy, meshYZy, meshYZz);
    pZX = PhysicalPlane(meshZXx, 0*meshZXz, meshZXz);
    
    %Сохраняем индексы для массивов значений
    pXY_Ex_idx = pXY.addValues;
    pXY_Ez_idx = pXY.addValues;
    pYZ_Ex_idx = pYZ.addValues;
    pYZ_Ez_idx = pYZ.addValues;
    pZX_Ex_idx = pZX.addValues;
    pZX_Ez_idx = pZX.addValues;

    grid1 = FlightView(ax1, params, Nparticles);    
    grid1 = grid1.addPlane(pZX, ColorPlane.phys2guiTranslate([0; params.visualYRange(1); 0]), pZX_Ex_idx);
    colorbar1 = colorbar(ax1);
    colorbar1.Position = ax1.Position;
    colorbar1.Position(1) = Column1Width + Column2Width;
    colorbar1.Position(3) = 0.95 - colorbar1.Position(1);
    title(colorbar1, 'Er, V/m');
    
    grid2 = FlightView(ax2, params, Nparticles);    
    grid2 = grid2.addPlane(pZX, ColorPlane.phys2guiTranslate([0; params.visualYRange(1); 0]), pZX_Ez_idx);
    colorbar2 = colorbar(ax2);
    colorbar2.Position = ax2.Position;
    colorbar2.Position(1) = Column1Width + Column2Width;
    colorbar2.Position(3) = 0.95 - colorbar2.Position(1);
    title(colorbar2, 'Ez, V/m');
    
   
    gridside1 = FlightView(axside1, params, Nparticles);
    gridside1 = gridside1.addPlane(pXY, ColorPlane.phys2guiTranslate([0; 0; params.visualZRange(1)]), pXY_Ex_idx);
    
    gridside2 = FlightView(axside2, params, Nparticles);
    gridside2 = gridside2.addPlane(pXY, ColorPlane.phys2guiTranslate([0; 0; params.visualZRange(1)]), pXY_Ez_idx);
    gridside2 = gridside2.addPlane(pYZ, ColorPlane.phys2guiTranslate([params.visualXRange(1); 0; 0]), pYZ_Ez_idx);
    gridside2 = gridside2.addPlane(pZX, ColorPlane.phys2guiTranslate([0; params.visualYRange(1); 0]), pZX_Ez_idx);
    
    
    %Настраиваем вид и камеру
    set(ax1, 'XLim', params.visualZRange, 'YLim', params.visualXRange);
%     set(ax2, 'XLim', visualZRange, 'YLim', visualXRange);
    set(axside1, 'XLim', params.visualXRange, 'YLim', params.visualYRange);
    set(axside2, 'XLim', params.visualXRange, 'YLim', params.visualYRange, 'ZLim', params.visualZRange);
%     set(ax1, 'Units', 'centimeters');    
    

    %метка времени
    ttt = text(ax1, 0.1, 1, 'time', 'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

    %Моменты времени для кадров
%     midt = min(t) + (max(t) - min(t))/100;
%     tFrames = [linspace(min(t), midt, 1 + 1e+1) linspace(midt, max(t), 1 + 1e+2)];
    tFrames = linspace(min(t), max(t), 1 + 1e+3);
    
    vidObj = VideoWriter(params.fileName);
    open(vidObj);
    %Анимация
    ax1Controller = AxesController(ax1, 'ZX');
    ax2Controller = AxesController(ax2, 'ZX');
    axside1Controller = AxesController(axside1, 'YX');
    axside2Controller = AxesController(axside2, '3D');
    for n = 1 : numel(tFrames)        
        tFrame = tFrames(n);
        [cx, cy, cz] = deal(CX(tFrame), CY(tFrame), CZ(tFrame));
        Xcur = interp1(t, X', tFrame)';
        Ycur = interp1(t, Y', tFrame)';
        Zcur = interp1(t, Z', tFrame)';
        grid1.setParticles(Xcur, Ycur, Zcur);
        grid2.setParticles(Xcur, Ycur, Zcur);
        gridside1.setParticles(Xcur, Ycur, Zcur);
        gridside2.setParticles(Xcur, Ycur, Zcur);                        
     
        grid1.moveCamera(0, 0, cz);
        grid2.moveCamera(0, 0, cz);
        gridside1.moveCamera(0, 0, cz);
        gridside2.moveCamera(0, 0, cz);
        
        
        pXY.translateFromStart([0; 0; cz]);
        pYZ.translateFromStart([0; 0; cz]);
        pZX.translateFromStart([0; 0; cz]);
        
        [Exy, Bxy] = EMfieldSolution.getFieldXYZt_EB(pXY.vertices, tFrame);
        pXY.setValue(pXY_Ex_idx, Exy(1, :));
        pXY.setValue(pXY_Ez_idx, Exy(3, :));
        [Eyz, Byz] = EMfieldSolution.getFieldXYZt_EB(pYZ.vertices, tFrame);
        pYZ.setValue(pYZ_Ex_idx, Eyz(1, :));
        pYZ.setValue(pYZ_Ez_idx, Eyz(3, :));
        [Ezx, Bzx] = EMfieldSolution.getFieldXYZt_EB(pZX.vertices, tFrame);
        pZX.setValue(pZX_Ex_idx, Ezx(1, :));
        pZX.setValue(pZX_Ez_idx, Ezx(3, :));
        
        grid1.updateColorPlanes;
        grid2.updateColorPlanes;
        gridside1.updateColorPlanes;
        gridside2.updateColorPlanes;
                
        CameraRotate(axside2, tFrame - min(t), max(t) - min(t));
                
        ttt.String = sprintf("t = %5.2f ps\n", tFrame*1e+12);
        
        ax1Controller.apply;
        ax2Controller.apply;
        axside1Controller.apply;
        axside2Controller.apply;
             
%         pause(0.1);
        caxis(axside1, caxis(ax1));
        caxis(axside2, caxis(ax2));
        drawnow;
        currFrame = getframe(fig);
        writeVideo(vidObj, currFrame);
    end
    
    close(vidObj);
end


function map = WongColorMap
%ColorMap как в статье Wong
    map = TripleColorMap([1; 0; 0], [1; 1; 1], [0; 0; 1]);
end

function map = MyColorMap
%ColorMap как в статье Wong
    map = TripleColorMap([0.8; 0.2; 0], [1; 1; 1], [0.2; 0; 0.8]);
end

function map = TripleColorMap(highColor, midColor, lowColor)
%модифицированный ColorMap
    N = 256;
    map = zeros(N, 3);
    for n = 1 : N
        level = 2*(n/N - 0.5);
        if level < 0
            level = -level;
            map(n, :) = lowColor*level + midColor*(1-level);
        else
            map(n, :) = highColor*level + midColor*(1-level);
        end
    end    
end


function CameraRotate(axside2, t, maxt)    
    phi = 10*pi*t/maxt;
    A = 0.2;
    r = [0; 1; 0] + [A; 0; 0]*cos(phi) + [0; 0; A]*sin(phi);    
    set(axside2, 'CameraPosition', r, 'CameraUpVector', [1; 0; 0]);
    box(axside2, 'on');
end
