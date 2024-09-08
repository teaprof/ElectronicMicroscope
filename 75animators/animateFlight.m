function animateFlight(EMfieldSolution, traj)
    %%Рисуем анимацию
    figure;
    hold on
    [X, Y, Z] = deal(traj.X, traj.Y, traj.Z);
    CX = mean(X, 1); CY = mean(Y, 1); CZ = mean(Z, 1);
    
    radius = EMfieldSolution.r(1);
    zmin = min(min(Z));
    zmax = max(max(Z));
    cyl = mycylinder([0; 0; zmin], [0; 0; zmax], radius);
    cyl.FaceAlpha = 0.2;
    
    xlim = [-radius, radius];
    ylim = [-radius, radius];
    zlim = [zmin, zmax];
%     zlim = [min(Z-CZ, [], 'all'), max(Z-CZ, [], 'all')];
%     minlim = min([xlim ylim zlim]);
%     maxlim = max([xlim ylim zlim]);
%     xlim = [minlim maxlim];
%     ylim = [minlim maxlim];
%     zlim = [minlim maxlim];
    set(gca, 'XLim', xlim, 'YLim', ylim, 'ZLim', zlim);
            
    sc = scatter3(r0(1, :), r0(2, :), r0(3, :), '.');
    grid('on');
    set(gca, 'CameraPosition', [xlim(2); ylim(2); zlim(2)]);
    tt = text(0, 0, 'ttt', 'Units', 'normalized');
    
    visualErPlane = mycolorplane([0 1 0], [0 0 1], xlim, zlim);
    
    t_avi = linspace(min(t), max(t), 100);
%     vidObj = VideoWriter(sprintf('booom_%g_MeV.avi', Estart*1e-6));
    vidObj = VideoWriter(sprintf('booom.avi'));
    open(vidObj);
    for k = 1 : numel(t_avi)
        XData = zeros(1, N);
        YData = zeros(1, N);
        ZData = zeros(1, N);
        for n = 1 : N
            XData(n) = interp1(t, X(n, :), t_avi(k));
            YData(n) = interp1(t, Y(n, :), t_avi(k));
            ZData(n) = interp1(t, Z(n, :), t_avi(k));
        end
        curCX = interp1(t, CX, t_avi(k));
        curCY = interp1(t, CY, t_avi(k));
        curCZ = interp1(t, CZ, t_avi(k));
        sc.XData = XData;
        sc.YData = YData;
        sc.ZData = ZData;
        tt.String = sprintf('%g ns', t_avi(k)*1e9);
%         pause(0.1);        
%         axis('equal');
        set(gca, 'XLim', xlim, 'YLim', ylim, 'ZLim', zlim);
        set(gca, 'CameraPosition', [4*radius; 0; curCZ]);
        set(gca, 'CameraTarget', [0; 0; curCZ]);
        set(gca, 'CameraViewAngle', 30);
        
        visualErPlane = WaveguidePlotEr(EMfieldSolution, visualErPlane, t_avi(k));
                
        drawnow;
        currFrame = getframe(gcf);        
        writeVideo(vidObj,currFrame);         
    end
    close(vidObj);
end

function visualEr = WaveguidePlotEr(sol, visualEr, t) %рисовалка        
    [N, K] = size(visualEr.r);
    for n = 1 : N
        for k = 1 : K
            XYZ = visualEr.r{N, K};
            [EE, HH] = sol.getFieldXYZt_EH(1, XYZ, t);
            visualEr.obj.CData(n, k) = norm(EE);            
        end
    end
end

