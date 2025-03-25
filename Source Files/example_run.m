% This test case tests the cdecoding time complexity

clear all;
addpath('Core functions','Tables');
tStart = tic;

M = 200;
Rate  = 1/2;

K = ceil(M*Rate); % infomation length 

maxIter = 10;

ldpc_info = LDPC_INFO(M,K); 


EsNodB = 0:0.25:0.25;
EsN0 = 10.^(EsNodB/10);
NoVec = 1./EsN0; % Es is normalized to 1
sigma = sqrt(NoVec/2);
[nBitErrs,nBlkErrs,BER_sim,FER_sim] = deal( zeros(size(EsNodB))); 

[audio, fs] = audioread('star1.wav');

audio_normalized = int16(audio * 32767);
audio_binary = dec2bin(typecast(audio_normalized(:), 'uint16'), 16);
binary_vector = audio_binary(:); 
binary_vector_int = binary_vector - '0';

t = 0;
for i = 1: length(EsNodB)% loop each SNR(EsN0)
    
    % print progress
    fprintf('\n Now running EsN0 %.1f dB [%d of %d]\n',EsNodB(i),i, length(EsNodB) );
    printLen = 0;
    tic;
    
    
    % loop each block
    Nblocks = 0;
    it = 1;
    decbinary = [ zeros(1,length(binary_vector))];

    while Nblocks < 7056  % stop criterion 
        
        %generate random A-bit message
        msg = binary_vector_int(it:it+K-1);
        %load('msg.mat'); msg = msg.';
        % append fillers to let it fit in base graph
        msg_ldpc = [msg; zeros(ldpc_info.n_F,1)];
        
        
        % encoding
        cword = nrldpc_encoder(ldpc_info,msg_ldpc);
        
        % rate matching
        txBits = nrldpc_rate_match(ldpc_info,cword,M);
        
        % BPSK bit to symbol mapping
        s = 1 - 2 * txBits; 
        %AWGN channel
        r = s + sigma(i) * randn(M,1); 
        %load('r.mat'); r = r.'; r = r(2*ldpc_info.Zc+1:end);
        % channel llr
      

        llr = 4 * EsN0(i) *r;
        %llr = r;
        % rate recovery
        llr = nrldpc_rate_recover(ldpc_info,llr);
        % decding 
        decBits = nrldpc_decoder(ldpc_info,llr,maxIter);
        decbinary(it:it+K-1) = decBits;
        it = it+K;

        %llr = llr(2*ldpc_info.Zc+1:end);
        %msg_cap = nrLDPCDecode(llr,ldpc_info.BGn,maxIter);%,'Algorithm','Offset min-sum');
        %decBits = msg_cap(1:K);
        %Counting errors
        Nerrs = sum(msg ~= decBits);

        if Nerrs > 0
            nBitErrs(i) = nBitErrs(i) + Nerrs;
            nBlkErrs(i) = nBlkErrs(i) + 1;
        end

        Nblocks = Nblocks + 1;
        BER_sim(i) = nBitErrs(i)/K/Nblocks;
        FER_sim(i) = nBlkErrs(i)/Nblocks;

        %print progress
        if mod(Nblocks,10) == 0  || Nblocks==1
            t = toc; 
            fprintf(repmat('\b',1,printLen));
            printStr = sprintf(' Elasped time is %.1f seconds,# Tx blocks: %d,# Error blocks: %d, BER: %.5f, BLER: %.5f',...
                                                                t,Nblocks,nBlkErrs(i),BER_sim(i),FER_sim(i));
            fprintf(printStr);
            printLen  = length(printStr);        
        end
    end


end
binaryStr = char(decbinary + '0');  % '1010'
binary_matrix = reshape(binaryStr, [], 16);
audio_integers = bin2dec(binary_matrix); 
audio_reconstructed = typecast(uint16(audio_integers), 'int16'); 
audio_reconstructed_normalized = double(audio_reconstructed) / 32767;



audiowrite('star2.wav', audio_reconstructed_normalized, fs);
scatterplot(s)
scatterplot(r)
    
% display and show simulation results    
fprintf('\n\n');
disp('simlulation results:');
disp('----------------------');
sim.M = ldpc_info.M;
sim.K = ldpc_info.K;
sim.R = Rate;

sim.EsNodB = EsNodB;
sim.BER = BER_sim;
sim.nBitErrs = nBitErrs;
sim.FER = FER_sim;
sim.nBlkErrs = nBlkErrs;
disp(sim);



