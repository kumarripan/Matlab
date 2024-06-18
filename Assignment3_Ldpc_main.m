clear all;
clc;

TbLen = 20496;
crc_poly = [1 1 0 0 0 0 1 1 0 0 1 0 0 1 1 0 0 1 1 1 1 1 0 1 1];% CRC polynomial for gCRC24A 38.212 section 5.1
L  = length(crc_poly) -1; % CRC Length
nPrb = 100;
nSym = 14;
nRE = 12;
nPilotRE = 6;
nREsPerPrb = (nSym* nRE) - nPilotRE; %NumOf dataREs per PRB
%N_totalRE = nPrb * nREsPerPrb;
modOrder = 2 ; % 2 bits for QPSK


TB = randi([0 1], 1, TbLen); %Generate ramdom bits of 0 &1 for TB length
crc_checksum = cal_crcChecksum(TB, crc_poly);%Calculate the CRC checksum
%disp('TB crc_checksum:');
%disp(crc_checksum);
TB_CRC = [TB, crc_checksum];% Append the generated  crc_checksum to TB "TB+CRC"

[Kcb, Kb, bg] = get_kcb_kb_bg(length(TB_CRC));
B = length(TB_CRC); % length of transport block + CRC
C = ceil(B/(Kcb - L)); % Number of codeBlocks (CB)
B_dash = B + (C*L);  % Total lentgh of all CB+CB_CRC
K_dash = B_dash/C;  % Size of each CB+CB_CRC
[Zc_value, iLS] = getZcValue(K_dash, Kb); % Get the ZC value and ILS index

G = nPrb * nREsPerPrb * modOrder;
E = G/C; %RateMatchedout per codeBlock
CB_CONCAT = 0;

if (Zc_value ~= 0)
    K  = Zc_value * Kb;
    F = K - K_dash; % Fillerbits

    CB_size = B/C;
    CB = zeros(C, CB_size);
    CB_CRC = zeros(C, CB_size + L);
    CB_CRC_Filler = zeros(C, CB_size + L + F);

    N = 0;
    % Tx side
    % 1: Generate the code Block = CB
    % 2: Append CRC to the CB => "CB+CRC"
    % 3: Append the Filler Bits to step-2 => "CB+CRC+Filler"
    % 4: Do LDPC encoding to each "CB+CRC+Filler"
    for i=1:C
        idx = i-1;
        startBit = idx * CB_size +1; %start bit index for CB 
        endBit = startBit + CB_size -1; % end bit index for CB
        CB(i,:) = TB_CRC(startBit: endBit); %code block
        crc = cal_crcChecksum(CB(i,:), crc_poly); %Generate CRC for CB              
        CB_CRC(i,:) = [CB(i,:), crc]; % Append CRC to CB => "CB+CRC"
        CB_CRC_Filler(i,:) = [CB_CRC(i,:), zeros(1, F)]; % Append Fillerbits =>"CB+CRC+Filler"

        %LDPC Encoding for each CB
        encodedBits(i,:) = LDPCEncode(transpose(CB_CRC_Filler(i,:)), bg); % LDPC encoding for "CB+CRC+Filler"
        N = length(encodedBits(i,:));
        rateMatchOut(i, :) = rateMatch(encodedBits(i,:), E);
        
        %codeBlockConcatinatination
        if  i == 1
            CB_CONCAT  =  rateMatchOut(i,:);
        else
            CB_CONCAT  = [CB_CONCAT, rateMatchOut(i,:)];
        end 
    end 

    %QPSK Modulation BLOCK    
    mod_output = qpskModulation(CB_CONCAT);
    %mod_output = 2*(CB_CONCAT-0.5); %BPSK Modulation
    scatterplot(mod_output);

    noise_power = (10^-5);
    noise = sqrt(noise_power)*randn(size(mod_output));
    rx_sig = mod_output + noise; % to add noise
    scatterplot(rx_sig);
        
    demod_out = qpskDemodulation(rx_sig);
    %demod_out = bpskDemodulation(rx_sig);

   % Receiver LDPC block 
    for i=1:C  
        CB_seg = demod_out(((i-1)*E)+1 : i*E);
        rateDematchOut = rateDematch(CB_seg, N);

        decodedBits_F(i,:) = LDPCDecode(transpose(rateDematchOut), bg, 25); 
        decodedBits_CRC(i,:) = [decodedBits_F(i, 1:K_dash)];

        error = numel(find(decodedBits_CRC(i,:) ~= CB_CRC(i,:)));
        disp("CB ERROR count");
        disp(error);

        crc_value = cal_crcChecksum(decodedBits_CRC(i,:), crc_poly);        
        if abs(crc_value) == 0
            disp('CB CRC PASSED');            
        else
            disp('CB CRC FAILED');
            disp(crc_value);
        end

        cbLength = length(decodedBits_CRC(i,:)) - L;
        if i ==1            
            TB_rx_crc = decodedBits_CRC(i,1:cbLength);
        else
            TB_rx_crc = [TB_rx_crc, decodedBits_CRC(i,1:cbLength)];
        end
    end
    TB_RX  = [TB_rx_crc(1: length(TB_rx_crc)-L)];
end


%disp('TB with CRC attached');    
%disp(TB_withCRC);

% Check the CRC is Passed or FAILED for the received TB+CRC
crc_value = cal_crcChecksum(TB_rx_crc, crc_poly);
if abs(crc_value) == 0
    disp('Received TB CRC PASSED');
else
    disp('Received TB CRC FAILED');
    disp(crc_value);
end
TB_error = numel(find(TB_RX ~= TB));
disp("Received TB ERROR count");
disp(TB_error);


% This function calculate the CRC checksum 
function checksum = cal_crcChecksum (TB, crc_poly)

    crc_poly_len = length(crc_poly);
    padding_len = crc_poly_len -1;
    padding_bits = zeros (1, padding_len);    
    
    A_padding = [TB, padding_bits]; % TB + Zero padding
    A_temp = A_padding;
    
        
    n_essentialbits = 0;
    essential_xor_out = 0;
    
    while (length(A_temp) >= crc_poly_len)
        
        xor_out = xor(A_temp(1:crc_poly_len), crc_poly);
        n_essentialbits = 0;
    
        for i = 1:crc_poly_len
            if 1 == (xor_out(1, i) & 1)
                essential_xor_out = xor_out(i:crc_poly_len);
                n_essentialbits = i;
                break;
            end
        end
    
        if n_essentialbits > 0
            A_temp = [essential_xor_out, A_temp(1, [(crc_poly_len +1):length(A_temp)])];
        else
            A_temp = A_temp(1, [(crc_poly_len +1):length(A_temp)]);
        end
        
    end
    
    if n_essentialbits == 0
        checksum = 0;
    else    
        checksum = A_temp;
        if (length(checksum) < (crc_poly_len -1))
            checksum = [zeros(1, crc_poly_len-length(checksum)-1), checksum];
        end
    end
end
    
%Fetch the ZC value
function [zcValue, iLS] = getZcValue(K_dash, Kb)
    % ZC table from 3gpp 38.212
    Zc= [2, 4, 8, 16, 32, 64, 128, 256;
     3, 6, 12, 24, 48, 96, 192, 384;
     5, 10, 20, 40, 80, 160, 320, 0;
     7, 14, 28, 56, 112, 224, 0, 0;
     9, 18, 36, 72, 144, 288, 0, 0;
     11, 22, 44, 88, 176, 352, 0, 0;
     13, 26, 52, 104, 208, 0, 0, 0;
     15, 30, 60, 120, 240, 0, 0, 0];

    zcValue = 0;
    iLS =  255;

    for i= 1:8
        for j = 1:8
            if ((Kb * Zc(i,j)) >= K_dash)
                zcValue = Zc(i,j);
                iLS = i -1;
                break;
            end
        end
    end
    
    if (zcValue == 0)
        disp("ERROR NO Zc Value found !!!")
    end
end
    
% Calculate the Kcb Kb and bg based on 3ggpp
% Kcb = maximum information block size = 22*384 for Bg1 & 10*384 for bg2
% Kb = number of systematic columns 22 for bg1 & 10 for bg2
% bg = Base Graph
function [Kcb, Kb, bg] = get_kcb_kb_bg(len)
    Kcb = 0;
    Kb = 0;
    bg = 0;

    % Below values from 3gpp standard
    KcbBg1 = 8448;
    KbBg1 = 22;
    KcbBg2 = 3840;
    KbBg2 = 10;     
    
    if (len >=KcbBg1)
        Kcb = KcbBg1;
        Kb = KbBg1;
        bg = 1;
    else
        Kcb  = KcbBg2;
        Kb = KbBg2;
        bg = 2;
    end
end
   
% RateMatching 
function rmOut = rateMatch(encodeBitsOut, E)
    L = length(encodeBitsOut);
    if L >= E
        rmOut = [encodeBitsOut(1:E)];
    else
        %e = L;
        %while e
        %    rmOut = [encodeBitsOut(1:lenEncodeOut)];
        %    e = E - L;
        %end
        disp (" NOT handled !!! rmOut < encodedBitsOut ")
    end

end

% Rate De-Matching 
function rdmOut = rateDematch(data, decodeLen)
    L = length(data);
    if decodeLen >= L
        rdmOut = data(1:L);
        rdmOut = [rdmOut, zeros(1,decodeLen-L)];
    else
       
        disp (" NOT handled !!! rdmOut < data length ")
    end

end

% Code Block Concatenation 
function codeBlockConcatinate(index, rateMatchOut, cbOut)
        if  index == 1
            cbOut  =  rateMatchOut;
        else
            cbOut  = [cbOut, rateMatchOut];
        end 
end

% QPSK Modulation
function modulation_data = qpskModulation(data)
    for i = 1:2:length(data)
        %I = 1/sqrt(2)* (1-2*data(i));
        %Q = 1/sqrt(2) * (1-2*data(i+1));
       
        I = 2*(data(i)-0.5);
        Q = 2*(data(i+1)-0.5);
        
        modOut = I + 1j*Q;
        if i == 1
            modulation_data = modOut;
        else
            modulation_data= [modulation_data, modOut];
        end
        %modOut = 2*(ldpc_coded_bits-0.5); %BPSK
    end
end

%QPSK demodulation
function demod_data = qpskDemodulation(data)
    for i= 1:length(data)
        I = real(data(i));
        Q = imag(data(i));

        llr_I0 = abs(-1 + I);
        llr_I1 = abs(1 + I);

        llr_I = log(llr_I0/llr_I1);

        llr_Q0 = abs(-1 + Q);
        llr_Q1 = abs(1 + Q);

        llr_Q = log(llr_Q0/llr_Q1);

        if i ==1
            %demod_data = [bit_I, bit_Q];
            demod_data = [llr_I, llr_Q];
        else
            %demod_data = [demod_data, bit_I, bit_Q];
            demod_data = [demod_data, llr_I, llr_Q];
        end
    end

end

%BPSK demodulation
function demod_data = bpskDemodulation(data)
        llr0 =  abs(-1 + data);   % BPSK demod
        llr1 =  abs(1 + data);    % BPSK demod
        
        llr = log(llr0./llr1);      % ldpc decoder requires log(p(r/0)/p(r/1))
        demod_data = llr;
end