classdef PhysicalPlane < handle
    properties
        meshx
        meshy
        meshz
        vertices0
        vertices
        values
    end
    methods
        function obj = PhysicalPlane(meshx, meshy, meshz)            
            assert(prod(size(meshx) == size(meshy)));
            assert(prod(size(meshx) == size(meshz)));
            obj.meshx = meshx;
            obj.meshy = meshy;            
            obj.meshz = meshz;
            
            Nrows = obj.rows;
            Ncols = obj.cols;            
            
            %create vertices
            obj.vertices = zeros(3, Nrows*Ncols);            
            for nrow = 1 : Nrows
                for ncol = 1 : Ncols
                    bottomleftidx  = obj.MultiIndexToSingleIndex(nrow, ncol);
                    obj.vertices0(:, bottomleftidx) = [meshx(nrow, ncol); meshy(nrow, ncol); meshz(nrow, ncol)];
                end
            end
            obj.vertices = obj.vertices0;
            obj.values = {};
        end
        function r = rows(obj)
            r = size(obj.meshx, 1);
        end
        function r = cols(obj)
            r = size(obj.meshx, 2);
        end
        function idx = MultiIndexToSingleIndex(obj, nrow, ncol)
            Nrows = obj.rows;
            idx = nrow + Nrows*(ncol-1);
        end
        function translateFromStart(obj, dr)
            obj.vertices = obj.vertices0 + dr;
        end
        function idx = addValues(obj)
            obj.values{end+1} = zeros(1, size(obj.vertices, 2));
            idx = numel(obj.values);
        end
        function setValue(obj, idx, values)
            obj.values{idx} = values;
        end
    end
end