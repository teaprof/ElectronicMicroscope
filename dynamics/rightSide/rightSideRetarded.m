function xt = rightSideRetarded(t, x, N, xhistory, thistory, EMfieldSolution)
% function xt = rightSideRetarded(t, x, N, xhistory, thistory, EMfieldSolution)
% The rhs for ode45 with retarded potential and interaction between particles 
%
% See also:
% rightSideNoInteraction, rightSideClassic, rightSideRetarded

%Несмотря на то, что учитывается запаздывающий потенциал, для уравнения
%движения используется нерелятивистское приближение, а именно: dv/dt = F.
%В релятивистком случае надо рассматривать d(gamma(t)*m*v(t))/dt = F,
%где gamma(t) = 1/sqrt(1 - beta^2), beta = v/c.
    q = getElectronCharge;
    m = getElectronMass;
    [r, v] = XtoPhase(x, N);
    rt = zeros(3, N);
    vt = zeros(3, N);
    for n = 1 : N
        rt(:, n) = v(:, n);
        nlayer = 2;
        [E, B] = EMfieldSolution.getFieldnXYZt_EB(nlayer, r(:, n), t);
        Scale = 1e+3;
        E = E * Scale;
        B = B * Scale*0;
        for k = 1 : N
            if n == k
                continue;
            end
            [xhist, vhist] = phaseToHistory(xhistory, k);
            [Ek, Bk] = getFieldRetarded(t, r(:, n),  xhist, vhist, thistory);
            E = E + Ek;
            B = B + Bk;
        end
        vt(:, n) = vt(:, n) + getAcceleration(q, m, v(:, n), E, B);
        assert(isreal(vt));
    end
    xt = [rt(:); vt(:)];
end

function [E, B] = getFieldRetarded(t, rdest, rsource_history, vsource_history, tsource_history)
    qsource = getElectronCharge;
    xsource_hist_fcn = griddedInterpolant(tsource_history, rsource_history(1, :));
    ysource_hist_fcn = griddedInterpolant(tsource_history, rsource_history(2, :));
    zsource_hist_fcn = griddedInterpolant(tsource_history, rsource_history(3, :));
    vxsource_hist_fcn = griddedInterpolant(tsource_history, vsource_history(1, :));
    vysource_hist_fcn = griddedInterpolant(tsource_history, vsource_history(2, :));
    vzsource_hist_fcn = griddedInterpolant(tsource_history, vsource_history(3, :));
    rsource_hist_fcn = @(t) [xsource_hist_fcn(t); ysource_hist_fcn(t); zsource_hist_fcn(t)];
    vsource_hist_fcn = @(t) [vxsource_hist_fcn(t); vysource_hist_fcn(t); vzsource_hist_fcn(t)];
    [E, B] = MovingChargeField.getEBanalytical(qsource, rsource_hist_fcn, vsource_hist_fcn, rdest, t);
end
