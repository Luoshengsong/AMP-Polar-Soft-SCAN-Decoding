classdef PolarSCAN
    properties
        design_snr_dB = 0

        FEC_Len % length of Polar coded length
        bitLen_beforeCoding
        
        FZlookup
        bitreversedindices
        F_kron_n
    end
    
    methods

        function obj = PolarSCAN(FEC_Len, bitLen_beforeCoding)
               obj.FEC_Len = FEC_Len;
               obj.bitLen_beforeCoding = bitLen_beforeCoding;
               [FZlookup1, bitreversedindices1, F_kron_n1] = init_Polar(obj);
               obj.FZlookup = FZlookup1;
               obj.bitreversedindices = bitreversedindices1;
               obj.F_kron_n = F_kron_n1;

           end
       function [FZlookup, bitreversedindices, F_kron_n] = init_Polar(obj)
           % N,K,n,construction_method,design_snr_dB,sigma,crc_size
            N = obj.FEC_Len;
            %K = obj.SourceCoding_Len + obj.CRC_Len;  % CRC bits included
            n = log2(N);
            F = [1 0;1 1];
            BB=1;
            for ii=1:n
                BB = kron(BB,F);
            end
            F_kron_n = BB;
            
            bitreversedindices = zeros(1,N);
            for index = 1 : N
                % for example ...
                % u_7 --> (7-1=6) --> (dec2->bin: 110) -->(flip:011) --> (bin->dec:3) --> (3+1=4) --> u_4
                bitreversedindices(index) = bin2dec(wrev(dec2bin(index-1,n)));
            end
            
            constructed_code_file_name = sprintf('constructedCode\\PolarCode_block_length_%d_designSNR_%.2fdB_method_BhattaBound.txt', N, obj.design_snr_dB);
            %should first use construct_polar_code(n) to construct the polar code
            indices = load(constructed_code_file_name);
            FZlookup = zeros(1,N);
            FZlookup(indices(1: obj.bitLen_beforeCoding)) = -1;
       end

       function Polar_Bits = Polar_SCAN_encode(obj, Appended_Bits)
            % (N,K) Polar codes
            % N = obj.FEC_Len; K = obj.SourceCoding_Len + obj.CRC_Len
            % Aims to map the CRC-Appended-Bits into Polar Bits
            [No_Active_Users, ~] = size(Appended_Bits);
            Polar_Bits = zeros(No_Active_Users, obj.FEC_Len);
            for active_user = 1: No_Active_Users
                x = obj.FZlookup;
                x (x == -1) = Appended_Bits(active_user, :); % -1's will get replaced by message bits below
                x = x(obj.bitreversedindices + 1);
                Polar_Bits(active_user, :) = mod(x * obj.F_kron_n, 2);
            end
       end

       function [u_llr, c_llr] = polar_SCAN_decode_alpha(obj, de_interlv_LLR_Bits, iter_num)
            % For (N, K) Polar code, we have
            N = obj.FEC_Len;
            K = obj.bitLen_beforeCoding;
            n = log2(N);
            
            Rate = K /N;
            if Rate <= 1/3
                alpha = 0.2;
            elseif Rate <= 1/2
                alpha = 0.6;
            elseif Rate <= 2/3
                alpha = 1.2;
            else
                alpha = 1.8;
            end
            
            plus_infinity = 1e3;
            Left_Msg = zeros(N,n+1); % left message
            Right_Msg = zeros(N,n+1); % Right message
            
            % initialization
            Left_Msg(:,n+1) = de_interlv_LLR_Bits';  % initial of left message
            Right_Msg(obj.FZlookup == 0, 1) = plus_infinity; % initial of right message
            
            for it = 1: iter_num
                for phi = 0:N-1
                    [Left_Msg, Right_Msg] = obj.updateLLRMap(n, phi, n, Left_Msg, Right_Msg);
                    if mod(phi,2)~=0
                        [Left_Msg, Right_Msg] = obj.updateBitMap(n, phi, n, Left_Msg, Right_Msg);
                    end
                end
            end
            
            mean_Right_Msg = mean(abs(Right_Msg(:,n+1)));
            mean_Left_Msg = mean(abs(Left_Msg(:,n+1)));
            %输出最终的左信息u_llr和右信息c_llr
            u_llr = Left_Msg(:,1)+Right_Msg(:,1);
            
            c_llr = Right_Msg(:,n+1)+alpha * mean_Right_Msg / mean_Left_Msg * Left_Msg(:,n+1);
            
            maxarg = 40;
            u_llr = max(min(u_llr,maxarg),-maxarg);
            c_llr = max(min(c_llr,maxarg),-maxarg);
        end
            
        function [Left_Msg,Right_Msg] = updateLLRMap(obj, lambda,phi,n,Left_Msg,Right_Msg)
            if lambda == 0
                return;
            end
            psi = floor(phi/2);
            if mod(phi,2)==0
                [Left_Msg,Right_Msg] = obj.updateLLRMap(lambda-1,psi,n,Left_Msg,Right_Msg);
            end
            for omega=0:2^(n-lambda)-1
                if mod(phi,2)==0
                    %do sth
                    Left_Msg( phi + omega * 2^lambda + 1, n + 1 - lambda ) = ...
                        obj.fFunction( ...
                        Left_Msg( psi + 2 * omega * 2^(lambda-1) + 1, n + 2 - lambda), ...
                        Left_Msg( psi + (2 * omega + 1) * 2^(lambda-1) + 1, n + 2 - lambda)+ ...
                        Right_Msg(phi+1+omega*2^lambda+1,n+1-lambda) ...
                        );
                else
                    %do sth
                    Left_Msg( phi + omega * 2^lambda + 1,n + 1 - lambda ) = ...
                        Left_Msg( psi + (2 * omega + 1) * 2^(lambda-1) + 1, n + 2 - lambda ) +...
                        obj.fFunction(...
                        Left_Msg( psi + 2 * omega * 2^(lambda-1) + 1, n + 2 - lambda ), ...
                        Right_Msg( phi - 1 + omega * 2^lambda + 1, n + 1- lambda )...
                        );
                end
            end
        end
            
            
        function [Left_Msg,Right_Msg] = updateBitMap(obj, lambda,phi,n,Left_Msg,Right_Msg)
        
            psi = floor(phi/2);
            if mod(phi,2)~=0
                for omega = 0:2^(n-lambda)-1
                    Right_Msg(psi+2*omega*2^(lambda-1)+1,n+2-lambda) = obj.fFunction(Right_Msg(phi-1+omega*2^(lambda)+1,n+1-lambda),Right_Msg(phi+omega*2^(lambda)+1,n+1-lambda)+Left_Msg(psi+(2*omega+1)*2^(lambda-1)+1,n+2-lambda));
                    Right_Msg(psi+(2*omega+1)*2^(lambda-1)+1,n+2-lambda) = Right_Msg(phi+omega*2^(lambda)+1,n+1-lambda) + obj.fFunction(Right_Msg(phi-1+omega*2^(lambda)+1,n+1-lambda),Left_Msg(psi+2*omega*2^(lambda-1)+1,n+2-lambda));
                end
                if mod(psi,2)~=0
                    [Left_Msg,Right_Msg] = obj.updateBitMap(lambda-1,psi,n,Left_Msg,Right_Msg);
                end
            end
        end
        
        function c = fFunction(obj, a, b)
            c = sign(a)*sign(b)*min(abs(a),abs(b));
        end
    end
end