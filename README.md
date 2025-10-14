1. Compile getEBClassicCoder_mex
2. addpath(genpath('.'))
3. optimMain


Use GPU
Use MATLAB interpreter
Use codegen: codegen getEBClassicCoder

benchEBClassicVariants

drawEM

simMain

optimMain

parpool before all

simMain approx 1 hour (2696 v3)

see results

see figures

gifsicle -O3 --color 256 -i input.gif -o output.gif

camelCase


rightSideClassic:
    %[E, B] = getFieldClassic(q, r, v);            %эталонная реализация
    %[E, B] = getFieldClassicVectorized(q, r, v);
    %[E, B] = getFieldClassicGPU(q, r, v);
    [E, B] = getFieldClassicCoderMex(q, r, v);     %самый быстрый, см. benchEBClassicVariants


retarded -- very slow


exportGraphics
