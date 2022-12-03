text = fileread('novel.txt'); % Read the text file
ascii_val = double(text);     % convert to double
my_text = dec2bin(ascii_val)-'0';  % Convert to 7 bit element binary array

% convert it into a single element array
b = my_text';  
my_text_new = b(:)';

% Source coding - Lempel Ziv encoding
% h = LZcoding(my_text_new); % Not done

%producing the Test_0.txt
outfile = fopen('Test_0.txt','w');
fprintf(outfile,'%s\n' ,h);
fclose(outfile);

% Channel coding - Hamming encoder
encoded = zeros(length(h),11);
% message length = 7
% code word length = 11 (using 4 parity)
for i=1:length(h)
    encoded(i,1:11) = hamming_encoder(h(i,1:7)); 
end

e_matrix = encoded';  % converting to an array
e = e_matrix(:)';

% writing the output as text_01
outfile = fopen('text_1.txt','w');
fprintf(outfile,'%i\n' , e);
fclose(outfile);

% Modulating the signal - 16 QAM
M = 16;
m = QAM_modulator( e, M, true);

% Transmission through channel
SNR = 10;
r = awgn (m, SNR , 'measured', 09);

% Digital demodulation
d = DemodQAM (r, M);

% Channel decoder
x = hammingcode_decoder(d);

% Source decoder
% t = LZdecoding(x);

% Output transducer - Converting back to  text
% dividing t array into 7-bit segments and converting to decimal
t = [];
r_ascii = [];
ascii = 0;

for n = 1:7:length(t)- 6
     for p = 0:1:6
        ascii = ascii + t(1, n+p)*2^(6-p); 
     end
     r_ascii(1, end+1) = ascii;
     ascii = 0;
end

%converting data into characters
RXtext = char(r_ascii);

%producing the Text_2.txt
outfile = fopen('Text_2.txt','w');
fprintf(outfile,'%s', RXtext);
fclose(outfile);

%Lempel Ziv decoder

% Hamming decoder
function decodedMessage = hammingcode_decoder( encodedMessage )

    n = length(encodedMessage);
    k = floor(log2(n));

    %populate a matrix of parity indices
    for i = 0:k
        k_arr(i+1) = 2^i; 
    end

    decodedMessage = [];
    checkBits = [];
    
    %this constructs a checkbit array to compare with the parity bits
    for i = 1:length(k_arr)
        parityBits(i) = encodedMessage(k_arr(i)); 
        paritySummer = 0;
        %begin at the parity bit index + 1
        j = k_arr(i)+1;
        %start counter at 1 to ensure parity bit isn't included in
        %calculation
        counter = 1;
        
        while (true)
            if counter >= k_arr(i)
                % skip k bits if k bits have been checked
                j = j + k_arr(i);
                %reset counter
                counter = 0;
            end
            
            if (j > length(encodedMessage))
                break;
            end
            
            paritySummer = paritySummer + encodedMessage(j);
            counter = counter + 1;
            j = j + 1;            
        end        
        parityBit = mod(paritySummer,2);
        checkBits(i) = parityBit;
    end
    
    %check for errors
    errBits = zeros(1,length(checkBits));
    errorFlag = 0;
    erroneousBit = 0;
    
    for i = 1:length(checkBits)
        if (checkBits(i) ~= parityBits(i))
           % record bits in error
           errBits(i) = 1;
           % set flag to 1 to indicate there were errors
           errorFlag = 1;
        end
    end
    
    if errorFlag == 0
        fprintf('No errors found!\n'); 
    else
        errors = find(errBits==1);
        %identify the bit in error
        erroneousBit = sum(k_arr(errors));        
    end
    
    % Construct error correct decoded message
    for i = 1:length(encodedMessage)
       if (~isempty(k_arr)) && i == k_arr(1)
           k_arr(1) = [];
       elseif i == erroneousBit
           %flip the erroneous bit and append to decoded message
           decodedMessage = [decodedMessage, ~erroneousBit];
       else
           %append data bit to decoded message
           decodedMessage = [decodedMessage, encodedMessage(i)];
       end
    end
    end

% Demodulating the signal - 16 QAM
function demodulated_signal = DemodQAM (RXsig, M)

    % Defining a zero array for complex numbers
    Dsig = complex(zeros(1,length(RXsig)));
    
    KMOD = 1/sqrt(10); % KMOD value for 64-QAM is 1/SQRT(42).
    RXsig = RXsig/KMOD; % Converting back to the original mapped values
    
    N = sqrt(M)-1; % 3
    % Checking each range to assign points to the points in constellation
    % Eg: points between 4<x<2 = 3
    
    for i = 1:length(RXsig)
        
        % real part
        found = false;
        for x = -N+1: 2: N-1  % check boundary values -6 to 6
            if real(RXsig(i)) < x
                Dsig(i) = x-1;
                found = true;
                break;
            end  
        end
       
        if ~found
           Dsig(i) = N; % check >6 values
        end
        
        % imaginary part
        found = false;
        for x = -N+1:2:N-1
            if imag(RXsig(i)) < x
                Dsig(i) = Dsig(i) + j*(x-1);
                found = true;
                break;
            end  
        end        
        if ~found
           Dsig(i) = Dsig(i) + j*N;
        end      
    end
    
    % Convert it into a bit stream
    % Defining the mapping table
    mapTable(1:4) = -3;     mapTable(5:8) = -1;
    mapTable(9:12) = +3;    mapTable(13:16) = +1;
    
    for i = 0:15
    if mod(i,4) == 0
      mapTable(i+1) = mapTable(i+1) -3*j;
    elseif mod(i+3,4) == 0
      mapTable(i+1) = mapTable(i+1) -1*j;
    elseif mod(i+1,4) == 0
      mapTable(i+1) = mapTable(i+1) +1*j;
    elseif mod(i+2,4) == 0
      mapTable(i+1) = mapTable(i+1) +3*j;
    end
    end
    
    % Defining a zero array for decimal index in the mapping table
     Dsig_index = zeros(1,length(Dsig));
     
     % Assigning index values by checking each mapping table points
     for x=1:length(Dsig)
         for y=1:length(mapTable)
             if Dsig(x) == mapTable(y) 
             Dsig_index(x) = y-1; 
             end
         end
     end  
    
    % converting to a binary array
    d_mat = (dec2bin(Dsig_index,4)- '0'); 
    d = d_mat';
    demodulated_signal = d(:)';
    
end

% Modulating the signal - 16 QAM
function Modulated_message = QAM_modulator( encodedMessage, M, avgPower)
    
    k = log2(M); % to get the number of bits in a binary symbol
    num = length(encodedMessage)/k; % Number of symbols

    mapTable(1:4) = -3;    mapTable(5:8) = -1;
    mapTable(9:12) = +3;   mapTable(13:16) = +1;
    
    for i = 0:15
    if mod(i,4) == 0
      mapTable(i+1) = mapTable(i+1) -3*j;
    elseif mod(i+3,4) == 0
      mapTable(i+1) = mapTable(i+1) -1*j;
    elseif mod(i+1,4) == 0
      mapTable(i+1) = mapTable(i+1) +1*j;
    elseif mod(i+2,4) == 0
      mapTable(i+1) = mapTable(i+1) +3*j;
    end
    end
    
    symbols = zeros(1, num); % Zero array for mappedSymbols
    
    for i = 1:k:length(encodedMessage)
        bin_symbol = encodedMessage(i:i+k-1);
        dec_value = 2^3*bin_symbol(1) + 2^2 * bin_symbol(2)...
                    + 2^1*bin_symbol(3) + 2^0*bin_symbol(4);
    
        symbols((i-1)/k+1) = mapTable(dec_value+1); 
    end
    
    %Construction of the signal
    m1 = symbols; 
    
    % Modify the function to have a Unit power output signal    
    KMOD = 1/sqrt(10); % KMOD value for 64-QAM is 1/SQRT(10).
    
    if avgPower == true
       Modulated_message = m1*(KMOD); 
    else
       Modulated_message = m1;       
    end


end

% Hamming encoder
function encodedMessage = hamming_encoder( binaryMessage )

    n = length(binaryMessage);

    % find highest binary multiplier
    i=0;
    while (true)
       if 2^(i) > n+(i)
           break;
       else
           i = i + 1;
       end
    end
    %highest binary power to identify largest index for parity bit
    k = i; 

    %create array of indices to inject parity bits
    for i = 0:k-1
       k_arr(i+1) = 2^i; 
    end

    encodedStrLen = n + length(k_arr);
    encodedMessage = 2*ones(1,encodedStrLen);

    for i = encodedStrLen:-1:1

        if (i == k_arr(end))
            paritySummer = 0;
            counter = 1;
            % start at first bit after parity bit
            j = i+1;
            while (true)
                % check to see if k bits have been checked
                if counter >= k_arr(end)
                    % skip k bits if k bits have been checked
                    j = j + k_arr(end);
                    %reset counter
                    counter = 0;
                end
                %check if at end of message
                if (j > length(encodedMessage))
                    break;
                end
                %cumulatively add the checked bits
                paritySummer = paritySummer + encodedMessage(j);
                counter = counter + 1;
                j = j + 1;
            end
            %take mod 2 of sum of checked bits
            parityBit = mod(paritySummer,2);
            encodedMessage(i) = parityBit;
            k_arr(end) = [];
        else
            %place message bit in new encoded message
            encodedMessage(i) = binaryMessage(end);
            binaryMessage(end) = [];
        end

    end
end

% Lempel ziv coder
