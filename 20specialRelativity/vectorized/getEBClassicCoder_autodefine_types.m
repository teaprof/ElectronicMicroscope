function getEBClassicCoder_autodefine_types
%function getEBClassicVectorized_autodefine_types
%Функция вызывает getEBClassicVectorized, позволяя Coder'у автоматически
%определить типы переменных-параметров
    N = 1e+4;
    qsource = getElectronCharge;
    r = rand(3, N);
    v = rand(3, N)*getSpeedOfLight;
    
    [E2, B2] = getEBClassicCoder(qsource, r, v);
    [E2, B2] = getEBClassicCoder(qsource, r(:, 1:N-1), v(:, 1:N-1));
end