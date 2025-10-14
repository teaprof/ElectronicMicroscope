function xt = rightSideNoInteraction(t, x, N, EMfieldSolution)
%function xt = rightSideNoInteraction(t, x, N, EMfieldSolution)
% The rhs for ode45 which doesn't takes into account the interaction between
% particles.
%
% See also:
% rightSideNoInteraction, rightSideClassic, rightSideRetarded
    q = getElectronCharge;
    m = getElectronMass;
    [r, v] = XtoPhase(x, N);
    [E, B] = EMfieldSolution.getFieldXYZt_EB(r, t);
    rt = v;
    vt = getAcceleration(q, m, v, E, B);
    assert(isreal(vt));
    xt = PhaseToX(rt, vt);
end
