function getEBClassicTest
% Test getEBClassicCoder, getEBClassicGPU, getEBClassicVectorized
% The ground-truth results are calculated by MovingChargeField.getEBclassic
    N = 1e+3;
    r = rand(3, N);
    v = rand(3, N);
    q = getElectronCharge;    

    % Prepare expected values
    expectedE = zeros(3, N);
    expectedB = zeros(3, N);
    for n = 1 : N
        [En, Bn] = MovingChargeField.getEBclassic(q, r, v, r(:, n));
        En(:, n) = 0;
        Bn(:, n) = 0;
        E(:, n) = sum(En, 2);
        B(:, n) = sum(Bn, 2);
        expectedE(:, n) = sum(En, 2);
        expectedB(:, n) = sum(Bn, 2);
    end
    %Test getEBClassicVectorized
    [receivedVectorizedE, receivedVectorizedB] = getEBClassicVectorized(q, r, v);        
    assert(max(abs(expectedE(:) - receivedVectorizedE(:))) < 1e-5)
    assert(max(abs(expectedB(:) - receivedVectorizedB(:))) < 1e-5)
    fprintf('[+] Vectorized code is OK\n')

    %Test getEBClassicCoder
    [receivedFromCoderE, receivedFromCoderB] = getEBClassicCoder(q, r, v);
    assert(max(abs(expectedE(:) - receivedFromCoderE(:))) < 1e-5)
    assert(max(abs(expectedB(:) - receivedFromCoderB(:))) < 1e-5)
    fprintf('[+] Code of mex-function is OK\n')

    %Test getEBClassicCoder_mex
    if exist('getEBClassicCoder_mex', 'file') == 3
        [receivedFromMexE, receivedFromMexB] = getEBClassicCoder_mex(q, r, v);   
        assert(max(abs(expectedE(:) - receivedFromMexE(:))) < 1e-5)
        assert(max(abs(expectedB(:) - receivedFromMexB(:))) < 1e-5)
        fprintf('[+] mex-function is OK\n')
    else
        fprintf('[ ] mex-function is not tested since it doesn''t exists\n')
    end

    %Test getEBClassicGPU
    if canUseGPU
        [receivedFromGPUE, receivedFromGPUB] = getEBClassicGPU(q, r, v);
        assert(max(abs(expectedE(:) - receivedFromGPUE(:))) < 1e-5)
        assert(max(abs(expectedB(:) - receivedFromGPUB(:))) < 1e-5)
        fprintf('[+] Code for GPU is OK\n')
    else
        fprintf('[ ] Code for GPU is not tested since GPU is unavailable\n')
    end
end