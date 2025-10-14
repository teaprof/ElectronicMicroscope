classdef WaveguideSolver
%classdef WaveguideSolverub
%   Solve the waveguide equation finding such values for kz for which the
%   solution exists. 
    properties
        sol
        solutions
    end
    methods
        function obj = WaveguideSolver(geom, w, m)
            obj.sol = WaveguideSolution(geom, w, m);
            obj.solutions = {};
        end
        
        function obj = solve(obj, nSolutions)
            if nargin < 2
                nSolutions = 1;
            end
            % Initialize feasible region [lb0, ub0]
            [x0, lb0, ub0] = obj.solToX;            
            lb0(lb0 == -Inf) = -1e+5;
            ub0(ub0 == Inf) = 1e+5;
            % Last elements of lb0 and ub0 corresponds to kz, we require
            % that the search range for kz in non-empty and will not be
            % collapsed (later we address to kz values as x(end, :):           
            % assert(lb0(end) < ub0(end))

            % Collapse the dimensions for which lb0 = ub0
            lb = WaveguideSolver.collapse(lb0, lb0, ub0);
            ub = WaveguideSolver.collapse(ub0, lb0, ub0);

            % Create target function
            f = @(x) targetFcn(obj, WaveguideSolver.uncollapse(x, lb0, ub0));

            % Set options for local optimizator algorithm
            opts = optimoptions('lsqnonlin', 'Display', 'none', 'FunctionTolerance', 1e-8, 'StepTolerance', 1e-12, 'MaxFunctionEvaluations', 1e+4, 'MaxIterations', 1e+4);

            % Select initial values for kz to run optimizer from
            kz_values = linspace(lb0(end), ub0(end), 1e+2);

            % Alternative options and kz values:
            %opts2 = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'none'); %'FunctionTolerance', 1e-8, 'StepTolerance', 1e-12, 'MaxFunctionEvaluations', 1e+6, 'MaxIterations', 1e+6);
            %kz_values = linspace(obj.sol.k_vacuum, 1, 1e+2);

            % Run optimizer starting from each kz
            xcur = zeros(numel(lb), numel(kz_values));
            ycur = zeros(1, numel(kz_values));
            parfor n = 1 : numel(kz_values)
                x0cur = obj.setKz(kz_values(n), x0, lb0, ub0);
                [xcur(:, n), ycur(n)] = lsqnonlin(f, WaveguideSolver.collapse(x0cur, lb0, ub0), lb, ub, opts);
                % For debug purposes: plot the solution
                %sol1 = obj.XtoSol(uncollapse(obj, xcur(:, n), lb0, ub0));
                %waveguidePlot(sol1);
            end

            % Remove bad solutions
            idx = ycur < 1e-5;
            xcur = xcur(:, idx);

            assert(sum(idx) > 0, "Can't find solution. Does it exist? Try to refine the set kz_values of initial values for kz")


            % Uncollapse
            x0cur = zeros(numel(lb0), size(xcur, 2));
            for n = 1 : size(xcur, 2)
                x0cur(:, n) = WaveguideSolver.uncollapse(xcur(:, n), lb0, ub0);                
            end

            % Join solutions that are close (we consider they equal)
            [~, idx] = WaveguideSolver.uniqapprox(x0cur(end, :), 1e-6, 1e-3);

            % Find best nSolutions solutions
            for n = 1 : min(nSolutions, numel(idx))
                obj.solutions{n} = obj.XtoSol(x0cur(:, idx(n)));
            end
            obj.sol = obj.solutions{1};
        end
                
        function y = targetFcn(obj, x)
            obj.sol = obj.XtoSol(x);
            % Use one of the following boundary conditions checkers:
            %y = getGUSdiscrepancySimple1(obj);
            %y = getGUSdiscrepancySimple2(obj);
            %y = getGUSdiscrepancySimple3(obj);
            y = getGUSdiscrepancy(obj);
            y = y*1e+10;
        end
        
        function [x, lb, ub] = setKz(obj, kz, x, lb, ub)
            x(end) = kz;
            lb(end) = kz;
            ub(end) = kz;
        end
        
        function [x, lb, ub] = solToX(obj)
            x = zeros(obj.sol.nlayers*2+1, 1);
            lb = zeros(obj.sol.nlayers*2+1, 1);
            ub = zeros(obj.sol.nlayers*2+1, 1);
            for n = 1 : obj.sol.nlayers
                x(n*2-1) = obj.sol.C1(n);
                x(n*2  ) = obj.sol.C2(n);
                
                lb(n*2-1) = obj.sol.C1span(n, 1);
                ub(n*2-1) = obj.sol.C1span(n, 2);
                lb(n*2  ) = obj.sol.C2span(n, 1);
                ub(n*2  ) = obj.sol.C2span(n, 2);
            end
            lb(end) = obj.sol.kzspan(1);
            ub(end) = obj.sol.kzspan(2);
            x(end) = obj.sol.kz;
        end
        
        function sol = XtoSol(obj, x)
            sol = obj.sol;
            for n = 1 : sol.nlayers                
                sol.C1(n) = x(n*2-1);
                sol.C2(n) = x(n*2);
            end
            sol.C1(1) = 1;
            sol.kz = x(end);
            sol = expand(sol);
        end

        
        %Граничные условия в общем виде записываются так:
        %1) на границе диэлектрик-диэлектрик:
        %       равенство тангенциальных составляющих полей E и H
        %2) на границе металл-диэлектрик: 
        %       равенство нулю тангенциальных составляющих E
        
       
        function delta = getGUSdiscrepancySimple1(obj)
        % Используется только функция getEzr
        % Hphi принимается пропорциональной D
            nlayers = obj.sol.nlayers;
            delta = zeros(nlayers*4 - 1, 1);
            for n = 1 : nlayers
                r = obj.sol.r(n);
                [Ezr1, dEzr1_dr] = getEzr(obj.sol, n, r);
                Hphi1 = obj.sol.eps(n)*dEzr1_dr;
                if n < nlayers
                    [Ezr2, dEzr2_dr] = getEzr(obj.sol, n+1, r);
                    Hphi2 = obj.sol.eps(n+1)*dEzr2_dr;
                else
                    %внешний слой - металл                    
                    Ezr2 = 0;
                    Hphi2 = Hphi1;
                end                
                %fprintf('n = %d\n', n);
                %fprintf('Ez2 - Ez1 = %d + 1i*%d\n', real(Ez2 - Ez1), imag(Ez2 - Ez1));
                %fprintf('Hphi2 - Hphi1 = %d + 1i*%d\n', real(Hphi2 - Hphi1), imag(Hphi2 - Hphi1));
                delta(n*4 - 3) = real(Ezr2 - Ezr1);
                delta(n*4 - 2) = imag(Ezr2 - Ezr1);
                delta(n*4 - 1) = real(Hphi1 - Hphi2);
                delta(n*4 - 0) = imag(Hphi1 - Hphi2);
            end
        end
        

        function delta = getGUSdiscrepancySimple3(obj)
            nlayers = obj.sol.nlayers;
            delta = zeros(nlayers*4 - 1, 1);
            phi = 0;
            for n = 1 : nlayers
                r = obj.sol.r(n);
                [Ezr1, Er1, Ephi1, Bz1, Br1, Bphiz1] = getComplexAmplitude_EH(obj.sol, n, r, 0);
                if n < nlayers
                    [Ezr2, Er2, Ephi2, Bz2, Br2, Bphiz2] = getComplexAmplitude_EH(obj.sol, n+1, r, 0);
                else
                    %внешний слой - металл                    
                    Ezr2 = 0;
                    Bphiz2 = Bphiz1;
                end                
                %fprintf('n = %d\n', n);
                %fprintf('Ez2 - Ez1 = %d + 1i*%d\n', real(Ez2 - Ez1), imag(Ez2 - Ez1));
                %fprintf('Hphi2 - Hphi1 = %d + 1i*%d\n', real(Hphi2 - Hphi1), imag(Hphi2 - Hphi1));
                delta(n*4 - 3) = real(Ezr2 - Ezr1);
                delta(n*4 - 2) = imag(Ezr2 - Ezr1);
                delta(n*4 - 1) = real(Bphiz2-Bphiz1);
                delta(n*4 - 0) = imag(Bphiz2-Bphiz1);
            end
        end

        function delta = getGUSdiscrepancy(obj)
            nlayers = obj.sol.nlayers;
            delta = zeros(nlayers*4 - 1, 1);
            phi = 0;
            for n = 1 : nlayers
                r = obj.sol.r(n);
                [Ez1, Er1, Ephi1, Hz1, Hr1, Hphi1] = getComplexAmplitude_EH(obj.sol, n, r, phi);
                %Граничные условия:
                %Для поля E: непрерывность касательных составляющих, а если
                %граница с металлом, то касательная составляющая равна
                %нулю.
                %Для поля B: непрерывность нормальных составляющих, если
                %граница с металлом, то нормальные составляющие равны нулю
                if n < nlayers
                    [Ez2, Er2, Ephi2, Hz2, Hr2, Hphi2] = getComplexAmplitude_EH(obj.sol, n+1, r, phi);
                else
                    %внешний слой - металл
                    Ez2 = 0; 
                    Ephi2 = 0;
                    Hz2 = Hz1;
                    Hphi2 = Hphi1;
                    %delta(n*4 - 3) = Ez2 - Ez1; 
                    %delta(n*4 - 2) = Ephi2 - Ephi1; %непрерывность касательных составляющих поля E, которое внутри металла равно 0
                    %delta(n*4 - 1) = Br2 - Br1; %магнитных зарядов не существует - удовл автоматом для монохроматичных полей        
                end                
                %fprintf('n = %d\n', n);
                %fprintf('Ez2 - Ez1 = %d + 1i*%d\n', real(Ez2 - Ez1), imag(Ez2 - Ez1));
                %fprintf('Hphi2 - Hphi1 = %d + 1i*%d\n', real(Hphi2 - Hphi1), imag(Hphi2 - Hphi1));
                delta(n*4 - 3) = real(Ez2 - Ez1);
                delta(n*4 - 2) = imag(Ez2 - Ez1);
                %delta(n*4 - 2) = real(Ephi2 - Ephi1);
                %delta(n*4 - 1) = real(Hz2 - Hz1);
                delta(n*4 - 1) = 1e+3*real(Hphi2 - Hphi1);
                delta(n*4 - 0) = 1e+3*imag(Hphi2 - Hphi1);
            end
        end
    end
    methods(Static)
        function [x, idx] = uniqapprox(x, abstol, reltol)
            [x, sortidx] = sort(x, 'descend');
            idx = 1;
            value = x(1);
            for n = 2 : numel(x)
                curvalue = x(n);
                if( abs(curvalue-value)/abs(value) < reltol || abs(curvalue-value)< abstol )
                    continue;
                end
                idx(end+1) = n;
                value = curvalue;
            end
            x = x(idx);
            idx = sortidx(idx);
        end        
        
        function xcollapsed = collapse(x, lb, ub)
            idx = lb ~= ub;
            xcollapsed = x(idx);
        end
        
        function xuncollapsed = uncollapse(x, lb, ub)
            xuncollapsed = lb;
            xk = 1;
            for k = 1 : numel(lb)
                if(lb(k) ~= ub(k))
                    xuncollapsed(k) = x(xk);
                    xk = xk + 1;
                end
            end
        end
    end    
end