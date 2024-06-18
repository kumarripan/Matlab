clear all;
clc;

%TB = [1 1 0 1 0 1 1 1]; % Transport Block
%crc_poly = [1 0 1]; % CRC polynomial for D^2 +1 

TbLen = 319784;
TB = randi([0 1], 1, TbLen); %Generate ramdom bits of 0 &1 for TB length
crc_poly = [1 1 0 0 0 0 1 1 0 0 1 0 0 1 1 0 0 1 1 1 1 1 0 1 1];% CRC polynomial for gCRC24A 38.212 section5.1

%disp('TB');
%disp(TB);

%Calculate the CRC checksum to be appended at end of TB
crc_checksum = cal_crcChecksum(TB, crc_poly); 
disp('crc_checksum:');
disp(crc_checksum);

% Attched the generated  crc_checksum at the end of TB to be transmitted
TB_withCRC = [TB, crc_checksum];
%disp('TB with CRC attached');    
%disp(TB_withCRC);

% Check the CRC is Passed or FAILED for the received TB+CRC
crc_value = cal_crcChecksum(TB_withCRC, crc_poly);

if abs(crc_value) == 0
    disp('CRC PASSED');
else
    disp('CRC FAILED');
    disp(crc_value);
end


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
    
