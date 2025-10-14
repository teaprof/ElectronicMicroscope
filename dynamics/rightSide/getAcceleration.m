function acc = getAcceleration(q, m, v, E, B)
    f = getEMforce(q, v, E, B);
    acc = f./m;
end

function f = getEMforce(q, v, E, B)
    fE = E*q;
    fM = q*cross(v, B);
    f = fE + fM;
end
