classdef ColorPlane
    properties
        physicalPlane  %phsyical vertices when initialized
        valuesIdx
        
        axes
        guiVertices %gui vertices in axes 'axes'
        faces
        patchObj
        
        %This used to calculate transformedVertices before drawing
        phys2guiTransform
    end
    methods
        function obj = ColorPlane(axes, physicalPlane, phys2guiTransform, valuesIdx)
            obj.phys2guiTransform = phys2guiTransform;            
            obj.axes = axes;
            obj.physicalPlane = physicalPlane;
            obj.valuesIdx = valuesIdx;
                        
            Nrows = physicalPlane.rows;
            Ncols = physicalPlane.cols;
            
            %create faces
            obj.faces = zeros((Nrows-1)*(Ncols-1), 4);
            for nrow = 1 : Nrows - 1
                for ncol = 1 : Ncols - 1
                    bottomleftidx  = physicalPlane.MultiIndexToSingleIndex(nrow,   ncol);
                    bottomrightidx = physicalPlane.MultiIndexToSingleIndex(nrow+1, ncol);
                    topleftidx     = physicalPlane.MultiIndexToSingleIndex(nrow,   ncol+1);
                    toprightidx    = physicalPlane.MultiIndexToSingleIndex(nrow+1, ncol+1);
                    obj.faces(nrow + (ncol-1)*(Nrows-1), :) = [bottomleftidx bottomrightidx toprightidx topleftidx];
                end
            end
            
            %Actual drawing
            obj = obj.updatePosition;
            S.vertices = obj.guiVertices';
            S.faces = obj.faces;
            S.FaceColor = 'interp';
            S.CData  = rand(size(obj.guiVertices, 2), 1);
            S.EdgeColor = 'none';
            
            obj.patchObj = patch(axes, S);
        end
        function obj = updateColor(obj)
            obj.patchObj.CData = obj.physicalPlane.values{obj.valuesIdx};
%             obj.patchObj.CData = rand(size(obj.guiVertices, 2), 1);
        end
        function obj = updatePosition(obj)
            aff = affine3d(obj.phys2guiTransform);
            obj.guiVertices = aff.transformPointsForward(obj.physicalPlane.vertices')';
            if numel(obj.patchObj) > 0
                obj.patchObj.Vertices = obj.guiVertices';
            end
        end                
    end
    methods(Static)
        function affine3dmatrix = phys2guiTranslate(translate)
            rot = eye(3);
            affine3dmatrix = [rot' [0; 0; 0]; translate' 1];
        end
        
        function affine3dmatrix = phys2gui(backtranslate, physAxesStr, guiAxesStr, translate)
            physAxes = parseAxes(physAxesStr);
            guiAxes  = parseAxes(guiAxesStr);
            rot = guiAxes*inv(physAxes);
            translate = -rot*backtranslate + translate;
            affine3dmatrix = [rot' [0; 0; 0]; translate' 1];
        end
        function DCM = parseAxes(str)
            DCM = zeros(3, 3);
            for n = 1 : numel(str)
                switch(str(n))
                    case 'x'
                        dir = [1; 0; 0];
                    case 'X'
                        dir = [-1; 0; 0];
                    case 'y'
                        dir = [0; 1; 0];
                    case 'Y'
                        dir = [0; -1; 0];
                    case 'z'
                        dir = [0; 0; 1];
                    case 'Z'
                        dir = [0; 0; -1];
                    otherwise
                        error('unknwon axis');
                end
            end
            DCM(:, n) = dir;            
        end
    end
end

