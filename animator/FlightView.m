classdef FlightView
    properties
        colorPlanes
        particles
        params
        axes
    end
    methods
        function obj = FlightView(axes, params, NParticles)
            obj.params = params;
            obj.axes = axes;
            Zeros = zeros(NParticles, 1);
            obj.particles = scatter3(axes, Zeros, Zeros, Zeros, 3, 'k');
            set(axes, 'XLim', params.visualXRange, 'YLim', ...
                params.visualYRange, 'ZLim', params.visualZRange);
            xlabel(axes, 'x');
            ylabel(axes, 'y');
            zlabel(axes, 'z');
            set(axes, 'CameraUpVector', [1; 0; 0]);
            set(axes, 'CameraPosition', [0; -1; 0]);
            obj.colorPlanes = {};
        end
        function obj = addPlane(obj, varargin)
                obj.colorPlanes{end+1} = ColorPlane(obj.axes, varargin{:});
        end
        function obj = updateColorPlanes(obj)
            for n = 1 : numel(obj.colorPlanes)
                obj.colorPlanes{n}.updateColor;
                obj.colorPlanes{n}.updatePosition;
            end
        end
        function obj = setParticles(obj, X, Y, Z)
            obj.particles.XData = X;
            obj.particles.YData = Y;
            obj.particles.ZData = Z;
        end
        function obj = moveCamera(obj, cx, cy, cz)
            set(obj.axes, 'XLim', cx + obj.params.visualXRange, ...
                'YLim', cy + obj.params.visualYRange, 'ZLim', cz + obj.params.visualZRange);
        end
    end
    methods(Static)
        function [meshx, meshy, meshz] = createMesh(params, axstr)
            gridX = linspace(params.visualXRange(1), params.visualXRange(2), params.gridSize(1));
            gridY = linspace(params.visualYRange(1), params.visualYRange(2), params.gridSize(1));
            gridZ = linspace(params.visualZRange(1), params.visualZRange(2), params.gridSize(2));
            switch(axstr)
                case 'x'                    
                    [meshy, meshz] = meshgrid(gridY, gridZ);
                    meshx = 0*meshz;
                case 'y'
                    [meshz, meshx] = meshgrid(gridZ, gridX);
                    meshy = 0*meshx;
                case 'z'
                    [meshx, meshy] = meshgrid(gridX, gridY);
                    meshz = 0*meshx;
            end
        end
    end
end
