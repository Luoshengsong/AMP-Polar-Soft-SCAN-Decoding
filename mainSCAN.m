clc; clear all;
rng('default');
warning('off');

basePath = [fileparts(mfilename('fullpath')) filesep];
addpath([basePath '/Polar_SCAN_decoding'])        % BiGAMP code

poolobj = gcp('nocreate'); delete(poolobj);
% % p = parpool(16);

%% Setting 
CRCenable = false;
SNRdB = 4: 2: 12;%[0, 4, 8, 12, 16];
monteCarlos = [1e2, 1e2, 1e2, 1e2, 1e3, 1e3, 1e3, 1e4, 1e4];
para.constell = [1, -1];  % BPSK modulation
para.outerIt = 6;
para.AMPIter = 15;
para.SCANIt = 5;

para.Ka = 16; % number of active users
M = 32;  % number of antennas
pathloss = ones(para.Ka, 1);

K0 = 16;
K_crc = 4;
if ~CRCenable
    K_crc = 0;
end
N = 64;
% Rate = (K0 + K_crc) / N;
para.K0 = K0;
para.K_crc = K_crc;
para.CRCenable = CRCenable;

%% Preporcessing
% interleaver and de-interleaver related: for row vector
interleaver_IDs = randperm(N);
[~, de_interleaver_IDs] = sort(interleaver_IDs);
interlv_func = @(Bits)  Bits(:, interleaver_IDs);  % Bits could be a matrix with size (No_Active_Users x FEC_Len)
de_interlv_func = @(interlv_Bits) interlv_Bits(:, de_interleaver_IDs);
para.interlv_func = interlv_func;
para.de_interlv_func = de_interlv_func;

if CRCenable
    CRC = CRC(K_crc);
end
Polar_SCAN = PolarSCAN(N, K0 + K_crc);

%% executation
% initialization
para.Prob_constells = ones(N, para.Ka, length(para.constell)) ./ length(para.constell);

BER = nan(length(SNRdB), 1);
FER = nan(length(SNRdB), 1);
for snr_id = 1: length(SNRdB)
    sigma2 = N / 10^(SNRdB(snr_id) / 10);
    para.sigma2 = sigma2;
    disp(['SNR = ' , num2str( SNRdB(snr_id) ) , 'dB.[', num2str(snr_id) ' of ', num2str(length(SNRdB)), ']']);
    tic,
    BER_SNR = nan(monteCarlos(snr_id), 1);
    FER_SNR = nan(monteCarlos(snr_id), 1);

    for mt = 1: monteCarlos(snr_id)
        %Polar_SCAN = PolarSCAN(N, K0 + K_crc);

        % --------------------- Transmitter side --------------------------
        % info bits
        u0 = randi([0, 1], K0, para.Ka);
        % append CRC
        if CRCenable
            u = CRC.encode(u0); % (K0 + K_crc) x Ka
        else
            u = u0;
        end
        % Polar encoding
        encodedBits = Polar_SCAN.Polar_SCAN_encode(u'); % Ka x N
        % interleaving
        inlv_bits = interlv_func(encodedBits);
        % BPSK modulation
        X = 2 * inlv_bits - 1;

        %---------------- Channel effect (BPSK, only I-path) --------------
        H = sqrt(1) * ( randn(M, para.Ka) );
        beta = sqrt( sum(abs(H).^2, 1) );
        H_norm = H ./ beta;
        
        % ---------------------- Receiver side ----------------------------
        Y = H * X + sqrt(1 * sigma2) * ( randn(M, N) );
        % perform AMP detect: output the likelihood information
        Info_Bits_hat = AMP_SCAN_Iter(Y, H_norm, beta, para, Polar_SCAN);

        % BER
        bit_comp = bitxor(Info_Bits_hat, u');
        BER_SNR(mt) = sum(bit_comp(:)) / para.Ka / N;
        
        % FER
%         [~, idx_check] = ismember(Info_Bits_hat, u', 'rows');
%         idx_check = idx_check(idx_check > 0);
        checkFER = sum(bit_comp, 2);
        FER_SNR(mt) = sum(checkFER > 0)  / para.Ka;
       
    end
    BER(snr_id) = sum(BER_SNR) / monteCarlos(snr_id),
    FER(snr_id) = sum(FER_SNR) / monteCarlos(snr_id),
    toc,
    disp('--------------------------------------');
end
% % delete(p);

%  SNRdB = -2:9
%  BER = [0.0896972656250000	0.0803613281250000	0.0639355468750000 0.0474804687500000	0.0320693359375000	0.0197470703125000	0.0102285156250000	0.00458964843750000	0.00175224609375000  5.245117187500000e-04	0.000127246093750000	2.02148437500000e-05]
% FER = [0.821875000000000	0.736250000000000	0.604375000000000	0.452500000000000 0.313312500000000	0.194750000000000	0.101750000000000	0.0463562500000000	0.0176000000000000 0.00523125000000000	0.00123125000000000	0.000225000000000000]

figure; semilogy( SNRdB, BER, 'm-o' ,...
              SNRdB, FER, 'b-s');
xlabel('SNR(dB)');
grid on;
legend('BER', 'FER')