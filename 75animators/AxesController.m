classdef AxesController
    properties
        axes
        mode 
    end
    methods
        function obj = AxesController(axes, mode, label)
            obj.mode = mode;
            obj.axes = axes;
            obj.apply;
        end
        function obj = apply(obj)
            XLim = obj.axes.XLim;
            YLim = obj.axes.YLim;
            ZLim = obj.axes.ZLim;
            switch(obj.mode)
                case 'XY'
                    cameraPosition = [mean(XLim) mean(YLim) 1];
                    cameraUpVector = [0; 1; 0];
                case 'YX'
                    cameraPosition = [mean(XLim) mean(YLim) 1];
                    cameraUpVector = [-1; 0; 0];
                case 'YZ'
                    cameraPosition = [1 mean(YLim) mean(ZLim)];
                    cameraUpVector = [0; 0; 1];
                case 'ZY'
                    cameraPosition = [-1 mean(YLim) mean(ZLim)];
                    cameraUpVector = [0; 1; 0];
                case 'ZX'
                    cameraPosition = [mean(XLim) 1 mean(ZLim)];
                    cameraUpVector = [1; 0; 0];
                case 'XZ'
                    cameraPosition = [mean(XLim) 1 mean(ZLim)];
                    cameraUpVector = [0; 0; 1];
                case '3D'
                    cameraPosition = [1; 1; 1];
                    cameraUpVector = [1; 0; 0];
            end
            set(obj.axes, 'CameraPosition', cameraPosition, 'cameraUpVector', cameraUpVector);
%             axes.ZDir = 'reverse';
%             axes.XAxisLocation = 'top';
        end
    end
end