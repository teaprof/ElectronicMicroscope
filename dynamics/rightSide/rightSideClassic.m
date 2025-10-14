function xt = rightSideClassic(t, x, N, EMfieldSolution, electronsPerParticle)
%function xt = rightSideClassic(t, x, N, EMfieldSolution, electronsPerParticle)
% The rhs for ode45 which takes into account the interaction between
% particles.
%
% See also:
% rightSideNoInteraction, rightSideClassic, rightSideRetarded

    q = getElectronCharge*electronsPerParticle;
    m = getElectronMass*electronsPerParticle;
    [r, v] = XtoPhase(x, N);
    rt = zeros(3, N);
    vt = zeros(3, N);

    rt = v;
    
    %Учитываем взаимодействие с внешним полем
    [E, B] = EMfieldSolution.getFieldXYZt_EB(r, t);
    vt = vt + getAcceleration(q, m, v, E, B);

    %Учитываем взаимодействие каждый-с-каждым
    %[E, B] = getFieldClassic(q, r, v);            %эталонная реализация
    %[E, B] = getFieldClassicVectorized(q, r, v);
    %[E, B] = getFieldClassicGPU(q, r, v);
    [E, B] = getFieldClassicCoderMex(q, r, v);     %самый быстрый, см. benchEBClassicVariants
    vt = vt + getAcceleration(q, m, v, E, B);
    assert(isreal(vt));
    
    %Формируем результат
    xt = [rt(:); vt(:)];
end


function [E, B] = getFieldClassic(q, r, v)
%Расчёт силы со сторны заряда с индексом source на заряд dest
    nparticles = size(r, 2);
    E = zeros(3, 1);
    B = zeros(3, 1);
    for n = 1 : nparticles
        [En, Bn] = MovingChargeField.getEBclassic(q, r, v, r(:, n));
        En(:, n) = 0;
        Bn(:, n) = 0;
        E(:, n) = sum(En, 2);
        B(:, n) = sum(Bn, 2);
    end
    assert(prod(~isnan(E(:))));
    assert(prod(~isnan(B(:))));
end

function [E, B] = getFieldClassicCoderMex(q, r, v)
    [E, B] = getEBClassicCoder_mex(q, r, v);
end


function [E, B] = getFieldClassicVectorized(q, r, v)
    [E, B] = getEBClassicVectorized(q, r, v);
end

function [E, B] = getFieldClassicGPU(q, r, v)
    [E, B] = getEBClassicGPU(q, r, v);
end
