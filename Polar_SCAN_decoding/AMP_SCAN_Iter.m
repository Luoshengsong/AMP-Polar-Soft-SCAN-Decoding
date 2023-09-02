function Info_Bits_hat = AMP_SCAN_Iter(Y, H_norm, beta, para, Polar_SCAN)

constell = para.constell;
sigma2 = para.sigma2;
Ka = para.Ka;
outerIt = para.outerIt;
AMPIter = para.AMPIter;
SCANIt = para.SCANIt;
Prob_constells = para.Prob_constells;  % N x Ka x length(constell)

Info_Bits_hat = nan( Ka, para.K0 + para.K_crc  );

% beta is a row vector
for it = 1: outerIt
    % AMP
    [X_effect_hat, X_effect_var] = AMP_detect(Y, H_norm, beta, constell, Prob_constells, sigma2, Ka, AMPIter);
    % Update X
    Xhat = diag(1./beta) * X_effect_hat ;
    Xvar = diag(1./ (beta.^2) ) * X_effect_var;

    % LLR conversion: BPSK
    LLR = ( - abs( Xhat + 1).^2 + abs( Xhat - 1).^2 ) ./ (2 * Xvar);
    LLR =  max(min(LLR, 40), -40);  % Ka x N

    % de-interleaving:
    deinlv_LLR = para.de_interlv_func(LLR);

    % decoding
    LLR_ext = zeros(size(deinlv_LLR));
    for k = 1: Ka
        [u_llr, c_llr] =Polar_SCAN.polar_SCAN_decode_alpha(deinlv_LLR(k,:), SCANIt);
        LLR_ext(k, :) = c_llr';
        if it == outerIt
            mhat_llr = u_llr(Polar_SCAN.FZlookup == -1, 1)';
            Info_Bits_hat(k, :) = (mhat_llr < 0);
        end
    end

    % interleaving
    LLR_ext = para.interlv_func(LLR_ext);

    % convert LLR into soft-symbol message
    Prob_bit0 = exp(LLR_ext) ./ (1 + exp(LLR_ext));
    Prob_bit1 = 1 - Prob_bit0;
    Prob_constells(:,:,1) = Prob_bit1';
    Prob_constells(:,:,2) = Prob_bit0';
end

% Info_Bits_hat = Info_Bits_hat';

end