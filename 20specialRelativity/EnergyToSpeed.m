function v = EnergyToSpeed(E_eV)
%function v = EnergyToSpeed(E_eV)
%Вычисляят скорость электрона [м/с] по его энергии [эВ]
    Erest = getElectronRestEnegry/abs(getElectronCharge); %eV
    assert(E_eV >= Erest);
    v = sqrt(1-(Erest/E_eV)^2)*getSpeedOfLight;        
end
