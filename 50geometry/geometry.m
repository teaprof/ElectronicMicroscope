classdef geometry
    properties
        %Слой номер i имеет радиусы:
        %   r(i-1) - внутренний, если i > 1, иначе 0
        %   r(i) - внешний
        %   eps(i) - диэлектрическую проницаемость
        %   mju(i) - магнитную проницаемость

        nlayers % количество слоёв
        r       % внешние радиусы слоёв, м
        eps     % комплексная диэлектрическая проницаемость слоёв
        mju     % комплексная магнитная проницаемость слоёв
    end
    methods
        %Constructor
        function obj = geometry(geomtype)
        %function obj = geometry([geomtype])
        % Initializes obj with the specified geometry and return it.
        % If geomtype is not provided create empty geometry. 
        % If geomtype is provided it can be one of the following:
        %   none - return empty object
        %   simple - return solid cylinder 
        %   waveguide - return cylindrical waveguide
        % These special types of geometry are used in tests
            if nargin == 0
                geomtype = '';
            end
            obj.nlayers = 0;
            obj.eps = [];
            obj.mju = [];
            obj.r = [];
            if strcmp(geomtype, 'none')
                return;
            elseif strcmp(geomtype, 'simple')
                %Содержит 100 слоёв из вакуума, решение должно совпадать с
                %однослойным решением для вакуума.
                nLayers = 100;
                for r = linspace(0, 10, nLayers + 1)
                    if r == 0
                        continue;
                    end
                    obj = addLayer(obj, 1, 1, r*1e-4);
                end
            elseif strcmp(geomtype, 'waveguide')                
                nLayers = 1;
                for r = linspace(0, 10, nLayers + 1)
                    if r == 0
                        continue;
                    end
                    obj = addLayer(obj, 1, 1, r*1e-4);
                end
                for r = linspace(0, 10, nLayers + 1)
                    if r == 0
                        continue;
                    end
                    obj = addLayer(obj, 5.5, 1, 1e-3 + r*1e-4); %todo: make as parameters
                end
            else
                assert(sprintf('Unknown geomtype: %s', geomtype))
            end
        end
        function obj = addLayer(obj, eps, mju, r_outer)
            obj.nlayers = obj.nlayers + 1;
            obj.eps(end+1) = eps;
            obj.mju(end+1) = mju;
            obj.r(end+1) = r_outer;
        end
    end    
end