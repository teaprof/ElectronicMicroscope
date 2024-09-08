classdef EMsimulator
    properties
        geom %структура с описанием геометрии и материалов задачи
        sol %структура с решением, описывающая ЭМ-поле
        mode %мода колебаний
        w %круговая частота 
    end    
    methods
        function obj = EMsimulator(geomtype)
            if nargin == 0
                geomtype = generateGeometry();
            else
                geom = generateGeometry('simple');
            end     
            w = 1e+14;
        end
        
        function obj = solve(obj)
            obj.sol = createSolStruct(obj.geom, obj.w);
            f = @(x) targetFcn(obj.sol, x);
            [x0, lb, ub] = solToX(obj.sol);
            [x, y] = lsqnonlin(f, x0, lb, ub);
            assert(y < 1e-14); 
            obj.sol = XtoSol(obj.sol, x);
        end
        
        function [Ez, Er, Ephi, Hz, Hr, Hphi] = getField(obj, nlayer, r, phi, sol)
            [Er, Ephi, Hr, Hphi] = getTransverse(nlayer, r, phi, obj.sol);
            Ez = obj.sol.Ezr{nlayer}(r)*obj.sol.Ezphi{nlayer}(phi);
            Hz = obj.sol.Hzr{nlayer}(r)*obj.sol.Hzphi{nlayer}(phi);
        end
                
        
        

    end
    
end
