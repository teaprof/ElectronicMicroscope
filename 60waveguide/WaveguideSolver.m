classdef WaveguideSolver
    properties
        sol
        solutions
    end
    methods
        function obj = WaveguideSolver(geom, w, m)
            obj.sol = WaveguideSolution(geom, w, m);
            obj.solutions = {};
        end
        
        function obj = solve(obj)
            [x0, lb, ub] = obj.solToX;
            lb(lb == -Inf) = -1e+5;
            ub(ub == Inf) = 1e+5;
            f = @(x) targetFcn(obj, unfold(obj, x, lb, ub));
            lbfolded = fold(obj, lb, lb, ub);
            ubfolded = fold(obj, ub, lb, ub);
            opts = optimoptions('lsqnonlin', 'Display', 'none', 'FunctionTolerance', 1e-8, 'StepTolerance', 1e-12, 'MaxFunctionEvaluations', 1e+4, 'MaxIterations', 1e+4);
            %opts2 = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'none'); %'FunctionTolerance', 1e-8, 'StepTolerance', 1e-12, 'MaxFunctionEvaluations', 1e+6, 'MaxIterations', 1e+6);
            %kz_values = linspace(obj.sol.k_vacuum, 1, 1e+2);
            kz_values = logspace(log10(ub(end)), log10(lb(end)), 1e+2); %last element of ub and lb is kz
            xcur = zeros(numel(lbfolded), numel(kz_values));
            ycur = zeros(1, numel(kz_values));
            parfor n = 1 : numel(kz_values)
            %for n = 1 : numel(kz_values)
                x0cur = obj.setKz(kz_values(n), x0, lb, ub);
                [xcur(:, n), ycur(n)] = lsqnonlin(f, fold(obj, x0cur, lb, ub), lbfolded, ubfolded, opts);
%                 sol1 = obj.XtoSol(unfold(obj, xcur(:, n), lb, ub));
%                 waveguidePlot(sol1);
%                 f(xcur(:, n));
            end
            %удаляем значения, для которых не удалось достичь околонулевых
            %значений ycur (суммы квадратов)
            idx = ycur < 1e-5;
            xcur = xcur(:, idx);
            [~, idx] = obj.uniqapprox(xcur(end, :)); %массив всех возможных kz
            for n = 1 : min(10, numel(idx))
                x = xcur(:, idx(n));
                obj.solutions{n} = obj.XtoSol(unfold(obj, x, lb, ub));
            end
            obj.sol = obj.solutions{1};
        end
        
        function [x, idx] = uniqapprox(obj, x)            
            [x, sortidx] = sort(x, 'descend');
            idx = 1;
            value = x(1);
            for n = 2 : numel(x)
                curvalue = x(n);
                if( abs(curvalue-value)/abs(value) < 1e-3 || abs(curvalue-value)< 1e-6 )
                    continue;
                end
                idx(end+1) = n;
                value = curvalue;
            end
            x = x(idx);
            idx = sortidx(idx);
        end        
        
        function xfolded = fold(obj, x, lb, ub)
            idx = find(lb ~= ub);
            xfolded = x(idx);
        end
        
        function xunfolded = unfold(obj, x, lb, ub)
            xunfolded = lb;
            xk = 1;
            for k = 1 : numel(lb)
                if(lb(k) ~= ub(k))
                    xunfolded(k) = x(xk);
                    xk = xk + 1;
                end
            end
        end
        
        function y = targetFcn(obj, x)
            obj.sol = obj.XtoSol(x);
%             y = getGUSdiscrepancySimple1(obj);
%             y = getGUSdiscrepancySimple2(obj);
%             y = getGUSdiscrepancySimple3(obj);
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
        %Используется только функция getEzr
        %Hphi принимается пропорциональной D
            nlayers = obj.sol.nlayers;
            delta = zeros(nlayers*4 - 1, 1);
            phi = 0;
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
%                 fprintf('n = %d\n', n);
%                 fprintf('Ez2 - Ez1 = %d + 1i*%d\n', real(Ez2 - Ez1), imag(Ez2 - Ez1));
%                 fprintf('Hphi2 - Hphi1 = %d + 1i*%d\n', real(Hphi2 - Hphi1), imag(Hphi2 - Hphi1));
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
%                 fprintf('n = %d\n', n);
%                 fprintf('Ez2 - Ez1 = %d + 1i*%d\n', real(Ez2 - Ez1), imag(Ez2 - Ez1));
%                 fprintf('Hphi2 - Hphi1 = %d + 1i*%d\n', real(Hphi2 - Hphi1), imag(Hphi2 - Hphi1));
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
%                     delta(n*4 - 3) = Ez2 - Ez1; 
%                     delta(n*4 - 2) = Ephi2 - Ephi1; %непрерывность касательных составляющих поля E, которое внутри металла равно 0
%                     delta(n*4 - 1) = Br2 - Br1; %магнитных зарядов не существует - удовл автоматом для монохроматичных полей        
                end                
%                 fprintf('n = %d\n', n);
%                 fprintf('Ez2 - Ez1 = %d + 1i*%d\n', real(Ez2 - Ez1), imag(Ez2 - Ez1));
%                 fprintf('Hphi2 - Hphi1 = %d + 1i*%d\n', real(Hphi2 - Hphi1), imag(Hphi2 - Hphi1));
                delta(n*4 - 3) = real(Ez2 - Ez1);
                delta(n*4 - 2) = imag(Ez2 - Ez1);
%                 delta(n*4 - 2) = real(Ephi2 - Ephi1);
%                 delta(n*4 - 1) = real(Hz2 - Hz1);
                delta(n*4 - 1) = 1e+3*real(Hphi2 - Hphi1);
                delta(n*4 - 0) = 1e+3*imag(Hphi2 - Hphi1);
            end
        end
        
        function plot(obj)
            waveguidePlot(obj.sol);
        end
    end
end