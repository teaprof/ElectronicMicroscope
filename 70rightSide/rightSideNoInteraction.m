function xt = rightSideNoInteraction(t, x, N, EMfieldSolution)
    q = getElectronCharge;
    m = getElectronMass;
    [r, v] = XtoPhase(x, N);
    rt = zeros(3, N);
    vt = zeros(3, N);
    [E, B] = EMfieldSolution.getFieldXYZt_EB(r, t);
    rt = v;
    vt = vt + getAcceleration(q, m, v, E, B);
    assert(isreal(vt));
    xt = PhaseToX(rt, vt);
end
