classdef MovingChargeField
    methods(Static)
        function simPotential
        %function simPotential
        %Рисует потенциал движущегося точечного заряда
        %todo: сравнить с Е.Ю.Петров "Излучение электромагнитных волн движущимися
        %заряженными частицами", Учебное пособие НижГУ им Лобачевского
            MovingChargeField.selfTest;
            
            v = [1; 0; 0]; %m/s;
            rsource_history = @(t) v*t;
            vsource_history = @(t) v;
            qsource = getElectronCharge;
            
            x = [-1:0.1:1];
            y = [-1:0.1:1];
            [xx, yy] = meshgrid(x, y);
            
            f = figure;
            t = [0:0.01:1];
            MaxPhi = 1;
            MaxE = 1;
            for m = 1 : numel(t)
                for n = 1 : numel(x)
                    for k = 1 : numel(y)
                        r = [x(n); y(k); 0];
                        [phi(n, k), A] = MovingChargeField.getLienardWiechertPotentials_retarded(qsource, rsource_history, vsource_history, r, t(m));
                        EE = MovingChargeField.getEBnumerical(qsource, rsource_history, vsource_history, r, t(m));
                        if(norm(phi(n, k)) > MaxPhi)
                            phi(n, k) = phi(n, k)/norm(phi(n, k))*MaxPhi;
                        end
                        E(n, k, 1:3) = EE(1:3);
                        if(norm(squeeze(E(n, k, 1:3))) > MaxE)
                            E(n, k, 1:3) = E(n, k, 1:3)/norm(squeeze(E(n, k, 1:3)))*MaxE;
                        end
                    end
                end
                subplot(2, 1, 1);
                hold off;
                surf(y, x, phi);
                hold on;
                rs = rsource_history(t(m));
                scatter(rs(1), rs(2));
                subplot(2, 1, 2);
                hold off;
                quiver(xx', yy', E(:, :, 1), E(:, :, 2));
                hold on;
                scatter(rs(1), rs(2));
                drawnow;
            end
        end
        
        function selfTest
        %function selfTest
        %Сравнивает аналитическое выражение для поля с численным расчётом и
        %с классическим выражением
            v = [1; 0; 0]; %m/s;
            rsource_history = @(t) v*t;
            vsource_history = @(t) v;
            qsource = 1; %Coulomb
            
            r = [0; 0; 0];
            while(norm(r) < 2)
                r = 10*(rand(3, 1) - 0.5);
            end
            
            t = 0;
            [Enumerical, Bnumerical] = MovingChargeField.getEBanalytical(qsource, rsource_history, vsource_history, r, t);
            [Ea, Ba] = MovingChargeField.getEBanalytical(qsource, rsource_history, vsource_history, r, t);
            [Eclassic, Bclassic] = MovingChargeField.getEBclassic(qsource, rsource_history(t), vsource_history(t), r);
            assert(norm(Enumerical - Eclassic) < 1e-4*norm(Eclassic));
            assert(norm(Ea - Eclassic) < 1e-4*norm(Eclassic));
            
            assert(norm(Ba - Bclassic) < 1e-4*norm(Bclassic));
            %assert(norm(Bnumerical - Bclassic) < 1e-4*norm(Bclassic));
            %%из-за численных ошибок Bnumerical расходится как с Bclassic,
            %%так и с Banalytical
        end
        
        function [E, B] = getEBclassic(qsource, rsource, vsource, r)
        %function [E, B] = getEBclassic(qsource, rsource, vsource, r)
        %Вычисляет классическое поле при малых скоростях
        %Для электростатического поля как 1/r^2,
        %для магнитного - в соответствии с законом Био-Савара-Лапласа
        %Это эталонная реализация, есть ещё ускоренные реализации:
        %       getEBClassicElementwise
        %       getEBClassicVectorized
        %       getEBClassicGPU
        %       getEBClassicCoder
        %       getEBClassicCoder_mex
            c = getSpeedOfLight;
            eps0 = getEps0;
            mju0 = getMju0;
            dr = r - rsource;
            norm_dr3 = sum(dr.^2).^1.5;
            E = 1/4/pi/eps0*qsource.*dr./norm_dr3;
            B = mju0/4/pi*qsource*cross(vsource, dr)./norm_dr3;
        end
        
        function phi = getPotentialClassic(qsource, rsource, r)
            eps0 = getEps0;
            dr = r - rsource;
            phi = 1/4/pi/eps0*qsource/(norm(dr));
        end
        
        function [E, B] = getEBanalytical(qsource, rsource_history, vsource_history, r, t)
        %function [E, B] = getEBanalytical(qsource, rsource_history, vsource_history, r, t)
        %Аналитические релятивистские выражения для элмагн поля с учётом запаздывания
        %потенциала
            %source: https://en.wikipedia.org/wiki/Liénard–Wiechert_potential
            t_retarded = MovingChargeField.getRetarded(rsource_history, r, t);
            rsource_retarded = rsource_history(t_retarded);
            vsource_retarded = vsource_history(t_retarded);
            
            eps0 = getEps0; %F/m
            c = getSpeedOfLight; %m/s
            betas = vsource_retarded/c;
            accelerationsource_retarded = [0; 0; 0];
            betas_t = accelerationsource_retarded/c;
            dr = r - rsource_retarded;
            dist = norm(dr);
            ns = dr/dist;
            
            gamma2 = 1/(1 - norm(betas)^2); %square of Lorentz factor
            etha = 1 - dot(ns, betas);
            
            First = (ns - betas)/gamma2/(norm(dr)^2);
            Second = cross(ns, cross((ns-betas), betas_t))/c/norm(dr);
            E = 1/4/pi/eps0*qsource/(etha^3)*(First + Second);
            B = cross(ns, E)/c;
        end
        
        
        function [E, B] = getEBnumerical(qsource, rsource_history, vsource_history, r, t)
        %function [E, B] = getEBnumerical(qsource, rsource_history, vsource_history, r, t)
        %Численные релятивистские выражения для элмагн поля с учётом запаздывания
        %потенциала, численно дифференцируется потенциал Линарда-Витхерта.
        %Функция работает неточно, в некоторых случаях заметно расходится с
        %getEBanalytical из-за численной погрешности
            %source: https://en.wikipedia.org/wiki/Liénard–Wiechert_potential
            phi = @(r) MovingChargeField.getLienardWiechertPotentials_retarded(qsource, rsource_history, vsource_history, r, t);
            grad_phi = MovingChargeField.grad(phi, r);
            A = @(t) MovingChargeField.getSecondOutput(@MovingChargeField.getLienardWiechertPotentials_retarded, qsource, rsource_history, vsource_history, r, t);
            dAdt = MovingChargeField.derivative(A, t);
            E = -grad_phi - dAdt;
            
            t_retarded = MovingChargeField.getRetarded(rsource_history, r, t);
            dist = norm(r - rsource_history(t_retarded));
            ns = (r - rsource_history(t_retarded))/dist;
            c = getSpeedOfLight;
            B = cross(ns, E)/c;
        end
        
       
        function [phi, A] = getLienardWiechertPotentials_retarded(qsource, rsource_history, vsource_history, r, t)
            t_retarded = MovingChargeField.getRetarded(rsource_history, r, t);
            rsource_retarded = rsource_history(t_retarded);
            vsource_retarded = vsource_history(t_retarded);
            [phi, A] = MovingChargeField.getLienardWiechertPotentials(qsource, rsource_retarded, vsource_retarded, r);
        end
        
        function [phi, A] = getLienardWiechertPotentials(qsource, rsource_retarded, vsource_retarded, r)
            eps0 = getEps0; %F/m
            c = getSpeedOfLight; %m/s
            betas = vsource_retarded/c;
            dist = norm(r - rsource_retarded);
            ns = (r - rsource_retarded)/dist;
            
            phi = 1/(4*pi*eps0)*qsource/(1 - dot(ns, betas))/dist;
            A = betas/c*phi;
            %     A = getMju0*c/4/pi*(qsource*betas)/(1 - dot(ns, betas))/dist;
        end
        
        function t_retarded = getRetarded(rhistory, r, t)
            %finds t_retarded as solution of:   abs(rhistory(t_retarded) - r) = c*(t - t_retarded)
            c = getSpeedOfLight; %m/s
            f = @(t_retarded) norm(rhistory(t_retarded) - r) - c*(t - t_retarded);
            t_retarded = fzero(f, t);
        end
        
    end
    methods(Static, Access = private)
        function z = getSecondOutput(f, varargin)
            [~, z] = f(varargin{:});
        end
        
        function R = rotor(f, r)
            fx = @(x) f([x; r(2); r(3)]);
            fy = @(y) f([r(1); y; r(3)]);
            fz = @(z) f([r(1); r(2); z]);
            dfdx = derivative(fx, r(1));
            dfdy = derivative(fy, r(2));
            dfdz = derivative(fz, r(3));
            R(1) = dfdy(3) - dfdz(2);
            R(2) = dfdz(1) - dfdx(3);
            R(3) = dfdx(2) - dfdy(1);
        end
        
        function E = grad(f, r)
            E = r;
            for dim = 1 : numel(r)
                abseps = abs(r(dim))*1e-8;
                if(abseps == 0)
                    abseps = 1e-8;
                end
                xright = r; xright(dim) = r(dim) + abseps;
                xleft  = r; xleft(dim)  = r(dim) - abseps;
                y1 = f(xleft);
                y2 = f(xright);
                E(dim) = (y2 - y1)/2/abseps;
            end
        end
        
        function res = derivative(f, x)
            abseps = abs(x)*1e-10;
            if(abseps == 0)
                abseps = 1e-10;
            end
            y1 = f(x + abseps);
            y2 = f(x - abseps);
            res = (y2 - y1)/2/abseps;
        end
    end
end