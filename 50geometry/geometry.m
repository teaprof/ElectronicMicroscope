classdef geometry
    properties
        nlayers %количество слоёв
        r %внешние радиусы слоёв, м
        eps %комплексная диэлектрическая проницаемость слоёв
        mju %комплексная магнитная проницаемость слоёв
        %Слой номер i имеет радиусы:
        %   r(i-1) - внутренний, если i > 1, иначе 0
        %   r(i) - внешний
        %диэлектрическую и магнитную проницаемости eps(i) и mju(i)
    end
    methods
        function obj = geometry(geomtype)
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
%                 obj = addLayer(obj, 1, 1, 1e-3);
                nLayers = 100;
                for r = linspace(0, 10, nLayers + 1)
                    if r == 0
                        continue;
                    end
                    obj = addLayer(obj, 1, 1, r*1e-4);
                end
            else
%                 obj = addLayer(obj, 5.5, 1, 1e-3);
%                 obj = addLayer(obj, 1, 1, 2e-3);
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
                    obj = addLayer(obj, 5.5, 1, 1e-3 + r*1e-4);
                end
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