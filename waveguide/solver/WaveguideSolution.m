classdef WaveguideSolution
%classdef WaveguideSolution
%   Представляет решение внутри коаксиального волновода, состоящего из двух
%   слоёв диэлектрика, обёрнутых металлической оболочкой
%       - ось z направлена по оси волновода
%   Порядок решения
%       Размерность вектора напряжённости электрического поля E: В/м
%       Ищем E в виде:
%           Ez(r, phi, z, t) = Ez(r)*exp(i*(w*t + kz*z + m*phi))
%           Ez(r) = C1*besselj(m, kxy*r) + C2*bessely(m, kxy*r), действительнозначное решение
%           phase = w*t + kz*z + m*phi;
%
%       Hzphi = 0
%       Hzr = 0
    properties(SetAccess = public, GetAccess = public)
        % Параметры геометрии задачи
        nlayers %к оличество слоёв
        r % внешние радиусы слоёв, м, numel(r) = nlayers
        eps % комплексная диэлектрическая проницаемость слоёв
        mju % комплексная магнитная проницаемость слоёв

        w         % угловая частота колебаний вынуждающего поля
        k_vacuum  % значение k в вакууме
        k2        % квадрат волнового вектора для каждого слоя
        
        %Решение
        m % мода по phi, для которой производится поиск решения
                        
        %Неизвестные коэффициенты в решении:
        kz  % проекция волнового вектора на kz
        C1  % коэффициенты при функциях Бесселя для каждого слоя
        C2  % коэффициенты при функциях Неймана для каждого слоя        

        %Фаза, рад
        phase % фаза в общем фазовом множителе
    end
    
    properties(SetAccess = private, GetAccess = public)
        % Зависимые коэффициенты в решении (обновляются с помощью функции-члена expand)
        kxy % проекция волнового вектора на плоскость xy (перп оси волновода)
    end

    properties(SetAccess = private, GetAccess = public)
        
        % Ограничения:
        % Ограничения учитываются при поиске решений для C1 и C2 численными
        % методами. 
        % C1span(i, :) = [a, b] - диапазон для C1(i)
        % C2span(i, :) = [a, b] - диапазон для C2(i)
        % kzspan(i, :) = [a, b] - диапазон для kz(i)
        % Если a = b, то значение для соответствующего коэффициента C
        % фиксировано и равно a.
        % Если a = b = 0, то член с соответствующей функцией Бесселя или Неймана не
        % добавляется в решение (это важно для внутреннего слоя, так как при r = 0
        % значение функции Неймана бесконечно)
        C1span % диапазон поиска решений по C1
        C2span % диапазон поиска решений по C2
        kzspan % диапазон для kz (kz  у всех слоёв имеет одинаковое значение)
    end
    methods
        function obj = WaveguideSolution(geometry, w, m)
        %function obj = WaveguideSolution(geometry, w, m)    
        % Инициализирует все поля
        % Поля C1 и C2 инициализируются начальным приближением
            c = getSpeedOfLight;
            %копируем все поля из geometry в sol
            obj.nlayers = geometry.nlayers;
            obj.eps = geometry.eps; 
            obj.mju = geometry.mju;
            obj.r = geometry.r;
            
            obj.k_vacuum = w/c;
            obj.k2 = obj.k_vacuum^2*(geometry.eps.*geometry.mju); %квадрат волнового вектора для каждого из слоёв            
                            
            obj.w = w;
            obj.m = m;
                
                                    
            %Диапазон для поиска kz:
            obj.kzspan = [1, min(sqrt(obj.k2) - 0.1)];
            assert(obj.kzspan(1) < obj.kzspan(2))            
            obj.kz = obj.kzspan(2);
            for n = 1 : geometry.nlayers
                obj.C1span(n, :) = [-Inf; Inf];
                obj.C2span(n, :) = [-Inf; Inf];    
                obj.C1(n) = 1; %коэффициент перед функцией Бесселя
                obj.C2(n) = 0; %коэффициент перед функцией Неймана
            end
            % требуем, чтобы этот коэффициент был равен 1, так как он будет 
            % задавать амплитуду решения (если учитывать только граничные
            % условия, то система  уравнений будет однородной => беск много решений)
            obj.C1span(1, :) = [1; 1]; 
            % требуем, чтобы для внутреннего слоя C2 = 0, так как функция 
            % Неймана не может быть решением (равна беск при r = 0)
            obj.C2span(1, :) = [0; 0]; 
  
          
            %Проверяем, чтобы начальное приближение попало в границы поиска
            assert(obj.kz >= obj.kzspan(1) && obj.kz <= obj.kzspan(2));
            for n = 1 : geometry.nlayers
                assert(obj.C1(n) >= obj.C1span(n, 1) && obj.C1(n) <= obj.C1span(n, 2));
                assert(obj.C2(n) >= obj.C2span(n, 1) && obj.C2(n) <= obj.C2span(n, 2));
            end
            
            obj.phase = 0;
        end
        
        function obj = expand(obj)
            %Вычисляет зависимые коэффициенты решения через независимые
            %в данном случае вычисляет kxy по формуле
            %       kxy = sqrt(k*2 - kz^2),
            %       k*2 = k_vacuum^2 * eps * mju
            for n = 1 : obj.nlayers
                kxy2 = obj.k_vacuum^2*obj.eps(n)*obj.mju(n) - obj.kz^2;
                if ~isreal(kxy2)
                    assert(isreal(kxy2));
                end
                obj.kxy(n) = sqrt(kxy2);
            end
        end
        
        function nlayers = findlayers(obj, XYZ)            
            nlayers = zeros(1, size(XYZ, 2));
            radius = sqrt(XYZ(1, :).^2 + XYZ(2, :).^2);
            rprev = 0;
            for n = 1 : numel(obj.r)
                idx = radius >= rprev & radius < obj.r(n);
                nlayers(idx) = n;
                rprev = obj.r(n);
            end
            nlayers(nlayers==0) = numel(obj.r) + 1;
        end
        
        
        function [E, B] = getFieldXYZt_EB(obj, XYZ, t)
            npoints = size(XYZ, 2);
            layers = findlayers(obj, XYZ);
            idx = layers > numel(obj.r);
            E = zeros(3, npoints);
            B = zeros(3, npoints);
            E(:, idx) = repmat([0; 0; 0], 1, sum(idx));
            B(:, idx) = repmat([0; 0; 0], 1, sum(idx));
            idx = layers <= numel(obj.r);
            [E1, B1] = getFieldnXYZt_EB(obj, layers(idx), XYZ(:, idx), t);
            E(:, idx) = E1;
            B(:, idx) = B1;
        end

        function [E, H] = getFieldXYZt_EH(obj, XYZ, t)
            nlayer = 1;
            cur_r = sqrt(XYZ(1)^2 + XYZ(2)^2);
            while(nlayer <= numel(obj.r) && obj.r(nlayer) < cur_r)
                nlayer = nlayer + 1;
            end
            if nlayer > numel(obj.r)
                E = [0; 0; 0];
                H = [0; 0; 0];
            else
                [E, H] = getFieldnXYZt_EH(obj, nlayer, XYZ, t);
            end
        end
        
        
        function [E, B] = getFieldnXYZt_EB(obj, nlayer, XYZ, t)
            [Ecplx, Hcplx] = getComplexAmplitudeXYZ_EH(obj, nlayer, XYZ, t);
            E = real(Ecplx);
            H = real(Hcplx);
            curmju = obj.mju(nlayer)*getMju0;
            B = H.*curmju;
        end

        function [E, H] = getFieldnXYZt_EH(obj, nlayer, XYZ, t)
            [Ecplx, Hcplx] = getComplexAmplitudeXYZ_EH(obj, nlayer, XYZ, t);
            E = real(Ecplx);
            H = real(Hcplx);
        end

        function [E, H] = getComplexAmplitudeXYZ_EH(obj, nlayer, XYZ, t)
        %Учитывается зависимость от r, phi и z.
            rr = sqrt(XYZ(1, :).^2 + XYZ(2, :).^2);
            phi = atan2(XYZ(2, :), XYZ(1, :));            
            [Ez, Er, Ephi, Hz, Hr, Hphi] = getComplexAmplitude_EH(obj, nlayer, rr, phi);
            Zexp = getPhaseMultiplierZT(obj, nlayer, XYZ(3, :), t);
            Ez = Ez.*Zexp; Er = Er.*Zexp; Ephi = Ephi.*Zexp;
            Hz = Hz.*Zexp; Hr = Hr.*Zexp; Hphi = Hphi.*Zexp;
            cosphi = cos(phi);
            sinphi = sin(phi);
            Ex = Er.*cosphi - Ephi.*sinphi;
            Ey = Er.*sinphi + Ephi.*cosphi;
            Hx = Hr.*cosphi - Hphi.*sinphi;
            Hy = Hr.*sinphi + Hphi.*cosphi;
            E = [Ex; Ey; Ez];
            H = [Hx; Hy; Hz]; %todo: всё-таки B или H?
        end

        function [Ez, Er, Ephi, Hz, Hr, Hphi] = getComplexAmplitude_EH(obj, nlayer, r, phi)
        %Если z-компоненты полей E и H имеют вид:
        %   Ez = Er(r, phi)*exp(i*(kz*z-w*t))
        %   Hz = Hr(r, phi)*exp(i*(kz*z-w*t))
        %то производная по z и по t в уравнения Максвелла может быть
        %выражена через умножение на i*kz или -i*w соотвественно, а
        %уравнения Максвелла с ротором можно рассматривать как СЛАУ и разрешить
        %относительно компонент по осям r и phi.
        %
        %В этой функции учитываются только зависимость компонент поля от r
        %и phi.
        %
        %Итоговые формулы взяты из Yeh Theory of Bragg fiber
        %
        %Ezr means Ez(r), Hzr means Hz(r)
        %Ezphi means Ez(phi), Hzphi means Hzphi(phi)
        %
        %For example, complex amplitude of z-component of field E:
        %        Ez = Ezr(r)*Ezphi(phi)*Ezz(z)

            
            curkz = obj.kz;
            cureps = obj.eps(nlayer)*getEps0;
            curmju = obj.mju(nlayer)*getMju0;
            curw = obj.w;
            
            [Ezr, dEzr_dr] = getEzr(obj, nlayer, r);
            [Ezphi, dEzphi_dphi] = getEzphi(obj, nlayer, phi);
            [Hzr, dHzr_dr] = getHzr(obj, nlayer, r);
            [Hzphi, dHzphi_dphi] = getHzphi(obj, nlayer, phi);
            
            Ez = Ezr.*Ezphi;
            Hz = Hzr.*Hzphi;

            % Commont coefficient in all formulas
            coef = 1i*curkz./(curw^2.*cureps.*curmju - curkz^2);
            
            % Calculate the partial derivatives for:
            %       Ez(r, phi) = Ezr(r)*Ezphi(phi)
            %       Hz(r, phi) = Hzr(r)*Hzphi(phi)            
            dEz_dr = dEzr_dr.*Ezphi;
            dEz_dphi = Ezr.*dEzphi_dphi;
            dHz_dr = dHzr_dr.*Hzphi;
            dHz_dphi = Hzr.*dHzphi_dphi;

            
            A = curw*curmju/curkz;
            %Calculate Er = coef*(dEz/dr + w*mju/kz*(1/r*dHz/dphi)):
            Er = dEz_dr + A./r.*dHz_dphi;
            Er = coef.*Er;
            Er(r == 0) = 0; %if r==0 then 0/0 = NaN occurs, but correct value is 0
            
            %Calculate Ephi = coef*(1/r*dEz/dphi - w*mju/kz*dHz/dr):
            Ephi = 1./r.*dEz_dphi - A.*dHz_dr;
            Ephi = coef.*Ephi;
            Ephi(r == 0) = 0; %if r==0 then 0/0 = NaN occurs, but correct value is 0
            
            %Calculate B:
            B = curw*cureps/curkz;
            %Calculate Hr = coef*(dHz/dr - w*eps/kz*(1/r*dEz/dphi)):
            Hr = dHz_dr - B./r.*dEz_dphi;
            Hr = coef.*Hr;
            Hr(r == 0) = 0; %if r==0 then 0/0 = NaN occurs, but correct value is 0
            
            %Calculate Hphi = coef*(1/r*dHz/dphi + w*eps/kz*dEz/dr):
            Hphi = 1./r.*dHz_dphi + B.*dEz_dr;
            Hphi = coef.*Hphi;
            Hphi(r == 0) = 0; %if r==0 then 0/0 = NaN occurs, but correct value is 0
        end
                
        function [Ezphi, dEzphi_dphi] = getEzphi(obj, nlayer, phi)
              Ezphi = cos(obj.m*phi); %если выбрать sin, то при m = 0 функция обратится в 0
              dEzphi_dphi = -sin(obj.m*Ezphi)*obj.m;
        end

        
        function [Ezr, dEzr_dr] = getEzr(obj, nlayer, r)
            Ezr = zeros(1, numel(r));
            dEzr_dr = zeros(1, numel(r));
            
            % Обрабатываем случай, если интервал поиска для C1 или C2 равен [0, 0]
            % Это означает, что соответствующая Бесселя или Неймана не должна фигурировать в решении
            BesselDisabled  = (obj.C1span(nlayer, 1) == 0 & obj.C1span(nlayer, 2) == 0);
            NeumannDisabled = (obj.C2span(nlayer, 1) == 0 & obj.C2span(nlayer, 2) == 0);
            BesselEnabled = ~BesselDisabled;
            NeumannEnabled = ~NeumannDisabled;
            
            BesselEnabledLayers = nlayer(BesselEnabled);
            Ezr(BesselEnabled) = Ezr(BesselEnabled) + obj.C1(BesselEnabledLayers).*besselj(obj.m, obj.kxy(BesselEnabledLayers).*r(BesselEnabled));
            dEzr_dr(BesselEnabled) = dEzr_dr(BesselEnabled) + obj.C1(BesselEnabledLayers).*besseljderiv(obj.m, obj.kxy(BesselEnabledLayers).*r(BesselEnabled)).*obj.kxy(BesselEnabledLayers);

            NeumannEnabledLayers = nlayer(NeumannEnabled);
            Ezr(NeumannEnabled) = Ezr(NeumannEnabled) + obj.C2(NeumannEnabledLayers).*bessely(obj.m, obj.kxy(NeumannEnabledLayers).*r(NeumannEnabled));
            dEzr_dr(NeumannEnabled) = dEzr_dr(NeumannEnabled) + obj.C2(NeumannEnabledLayers).*besselyderiv(obj.m, obj.kxy(NeumannEnabledLayers).*r(NeumannEnabled)).*obj.kxy(NeumannEnabledLayers);
        end
            
        
        function [Hzphi, dHzphi_dphi] = getHzphi(obj, nlayer, phi)
        %function [Hzphi, dHzphi_dphi] = getHzphi(obj, nlayer, phi)
        % Complex amplitude of z-component of magnetic field strength H:
        %           Hz = Hzr(r)*Hzphi(phi)*Hzz(z)
        % Since in the considered cases Hz is always zero it means that
        %              Hz(r) = 0, Hz(phi) = 0
            Hzphi = zeros(1, numel(phi));
            dHzphi_dphi = zeros(1, numel(phi));
        end        

        function [Hzr, dHzr_dr] = getHzr(obj, nlayer, r)
        %function [Hzphi, dHzphi_dphi] = getHzphi(obj, nlayer, phi)
        % Complex amplitude of z-component of magnetic field strength H:
        %           Hz = Hzr(r)*Hzphi(phi)*Hzz(z)
        % Since in the considered cases Hz is always zero it means that
        %              Hz(r) = 0, Hz(phi) = 0
            Hzr = zeros(1, numel(r));
            dHzr_dr = zeros(1, numel(r));
        end
        
        function phase = getPhaseMultiplierZT(obj, nlayer, z, t)
            phase = exp(1i*(obj.kz.*z - obj.w.*t + obj.phase));
        end
    end
end
