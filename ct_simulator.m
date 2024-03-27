%% ct-simulator
% The ct-simulator is a MATLAB script used to characterize the Packet Error
% Rate (PER) of radio transmissions using a Frequency-Shift Keying (FSK) 
% modulation scheme when two concurrent transmitters (ct) overlap in the 
% air sending the same bitstream in the presence of Additive White Gaussian 
% Noise (AWGN). It provides practical insights to understand the 
% performance of flooding protocols based on concurrent transmissions.
% 
% Author: Antonio Escobar <antonio@rednodelabs.com>

%% Background
% The main error source of concurrent transmissions is the beating effect. 
% When two transmissions of amplitudes _A1_ and _A2_ and frequencies _f1_ 
% and _f2_ overlap in the air, the superposition of both waves shows an 
% amplitude modulation (assuming they are perfectly synchronized at the 
% symbol level, i.e. triggered at the same time and neglecting path and 
% clock differences):
%
% $$ (A_1-A_2)+2A_2\left|\cos\left(2\pi\frac{f_{beat}}{2}t\right)\right| $$
% 
% Where _fbeat_ is the beating frequency:
%
% $$ f_{beat}=\left|f_1-f_2\right| $$

%% Comparing simulation results with well-known analytical figures
% For simplicity and easiness of comparison with analytical results, we
% consider noncoherent reception of binary FSK transmissions. 

% Modulation parameters
M = 2;         % Modulation order
Rs = 1e6;      % Symbol rate (1Mbps)
freqsep = Rs;  % Frequency separation (optimal for noncoherent detection)
nsamp = 8;     % Number of samples per symbol
Fs = Rs*nsamp; % Sample rate (Hz)

% Simulation parameters
EbNo = 0:2:12; % Eb/No (dB), desired simulation range 
N = 300000;    % Transmitted Bits
pktlen = 1;    % Packet length, in bits (with preamble)
preamblen = 0; % Preamble length, in bits

%%%
% First of all, we perform a basic sanity check, by setting to zero the
% amplitude of the second transmitter, effectively removing it from the
% simulation, and fixing the packet length to 1 bit. In that case, the
% result should be the well-known Bit Error Rate (BER) for noncoherent FSK
% receivers:
%
% $$ BER_{BFSK}=\frac{1}{2}\exp\left(-\frac{E_{b}}{2N_{0}}\right) $$

% Beating parameters
A2 = 0;        % Amplitude of transmitter 2, A1 is assumed to be 1
fbeat = 2e3;   % Beating frequency (Hz)

% Running the simulation
disp('Calculating BER for 1 Transmitter...')
ber_1 = ct_sim_run(M,Fs,nsamp,freqsep,A2,fbeat,EbNo,N,pktlen,preamblen);

% Analytical results for comparison
EbNo_linear = 10.^(EbNo/10);
ber_1_analytical = 0.5*exp(-EbNo_linear/2);

% Plotting the results
figure
plot(EbNo,ber_1,'b*')
hold on;
plot(EbNo,ber_1_analytical,'b-')
set(gca,'YScale','log')
xlabel('Eb/N0')
ylabel('BER')
legend('Simulated BFSK BER A2=0','Analytical BFSK BER A2=0')

%%%
% If both transmitters are received with the same energy (_A2_ = _A1_), the
% BER can also be obtained analytically*:
% 
% $$ BER=\frac{1}{2}\exp\left(-\frac{E_{b}}{N_{0}}\right)
% I_{0}\left(-\frac{E_{b}}{N_{0}}\right) $$
% 
% <https://publications.rwth-aachen.de/record/791687 *>

% Beating parameters
A2 = 1;        % Amplitude of transmitter 2, A1 is assumed to be 1
fbeat = 2e3;   % Beating frequency (Hz)

% Running the simulation
disp('Calculating BER for 2 Concurrent Transmitters...')
ber_2 = ct_sim_run(M,Fs,nsamp,freqsep,A2,fbeat,EbNo,N,pktlen,preamblen);

% Analytical results for comparison
ber_2_analytical = 0.5*exp(-EbNo_linear).*besseli(0,-EbNo_linear);

% Plotting the results
figure
plot(EbNo,ber_2,'b*')
hold on;
plot(EbNo,ber_2_analytical,'b-')
set(gca,'YScale','log')
xlabel('Eb/N0')
ylabel('BER')
legend('Simulated BFSK BER A2=1','Analytical BFSK BER A2=1')

%% Adding energy deltas
%%%
% In a real scenario, both ct will not be received with the same energy. If
% the relative amplitude of the second transmitter is reduced, i.e. the
% energy of one of the transmissions dominates the other on the receiver
% side, the impact of ct in the BER is reduced.

% Beating parameters
A2 = 0.5;      % Amplitude of transmitter 2, A1 is assumed to be 1
fbeat = 2e3;   % Beating frequency (Hz)

% Running the simulation
disp('Calculating BER for 2 Concurrent Transmitters with energy delta...')
ber_2_en_delta = ct_sim_run(M,Fs,nsamp,freqsep,A2,fbeat,EbNo,N,pktlen, ...
    preamblen);

% Plotting the results
figure
plot(EbNo,ber_2,'-*')
hold on;
plot(EbNo,ber_2_en_delta,'-*')
plot(EbNo,ber_1,'-*')
set(gca,'YScale','log')
xlabel('Eb/N0')
ylabel('BER')
legend('BFSK BER A2=1','BFSK BER A2=0.5','BFSK BER A2=0')

%% Analyzing packet errors
%%%
% If we analyze the error at the packet level, the beating frequency
% becomes relevant. The duration of the beating relative to the duration 
% of the packet determines how many destructive energy valleys a packet 
% transmission experiences, and how wide these valleys are. This has a 
% large impact in the Packet Error Rate (PER), since bit errors tend to 
% appear during the valleys*. We can analyze two different scenarios, 
% one in which the beating period is "wide" (longer than the packet 
% duration) and another in which it is "narrow" (shorter than the packet 
% duration).
%
% <https://dl.acm.org/doi/10.1145/3604430 *>
% 
% We aribitrarily choose a packet length of 128 bits and a beating period 
% twice as long as the packet duration to illustrate wide beating and a 
% beating period twice shorter for narrow beating. 

pktlen = 128;                    % Packet length, in bits (with preamble)
pduration = pktlen/(Rs*log2(M)); % Packet duration in the air, in seconds

% Beating parameters
A2 = 1;                % Amplitude of transmitter 2, A1 is assumed to be 1
fbeat = 0.5/pduration; % Beating frequency (Hz)

% Running the simulation
disp('Calculating PER for 2 Concurrent Transmitters with wide beating...')
[per_wide, errindx_wide] = ct_sim_run(M,Fs,nsamp,freqsep,A2,fbeat, ...
    EbNo,N,pktlen, preamblen);

% Beating parameters
A2 = 1;                % Amplitude of transmitter 2, A1 is assumed to be 1
fbeat = 2/pduration;   % Beating frequency (Hz)

% Running the simulation
disp(['Calculating PER for 2 Concurrent Transmitters ' ...
    'with narrow beating...'])
[per_narrow, errindx_narrow] = ct_sim_run(M,Fs,nsamp,freqsep,A2, ...
    fbeat,EbNo,N,pktlen, preamblen);

% Plotting the results
figure
plot(EbNo,per_wide,'-*')
hold on;
plot(EbNo,per_narrow,'-*')
set(gca,'YScale','log')
xlabel('Eb/N0')
ylabel('PER')
legend('BFSK PER (wide beating)','BFSK PER (narrow beating)')

%% Analyzing error indexes
%%%
% When a new packet reception starts, initial phase of the beating is 
% random, since it is determined by the relative frequency and phase 
% relationship between two independent oscillators. Therefore, a flat error 
% distribution is expected if we analyze the indexes of the errors across
% multiple packet receptions. Nevertheless, experimentally, if we plot the
% index of the errors after multiple packet receptions, we can see the 
% beating pattern*. This is due to the fact that for some packets (those in 
% which the preamble is broken) the reception will not even start. The 
% effect is a bias (beating) in the perceived error pattern, in which the 
% beating frequency can be observed.
%
% <https://dl.acm.org/doi/10.1145/3604430 *>
% 

% We choose an arbitrary preamble length of 32 bits (we assume errors in 
% the preamble will not be reported by the receiver)
preamblen=32;          % Preamble length, in bits

% Beating parameters
A2 = 1;                % Amplitude of transmitter 2, A1 is assumed to be 1
fbeat = 0.5/pduration; % Beating frequency (Hz)

% Running the simulation
disp(['Calculating PER for 2 Concurrent Transmitters with wide ' ...
    'beating and introducing the effect of the preamble...'])
[per_wide_preamb, errindx_wide_preamb,pktnorx_wide_preamb] = ...
    ct_sim_run(M,Fs,nsamp,freqsep,A2,fbeat,EbNo,N,pktlen, preamblen);

% Plotting the results
figure
EbNo_idx = ceil(length(EbNo/2)); %Plot for an arbitrary EbNo index
plot(1:length(errindx_wide(EbNo_idx,:)),errindx_wide(EbNo_idx,:),'-*')
hold on;
plot(1:length(errindx_wide_preamb(EbNo_idx,:)), ...
    errindx_wide_preamb(EbNo_idx,:),'-*')
xlabel('Bit Index')
ylabel('Error Count')
msg = sprintf(['wide beating (with preamble, %u%% of packets not ' ...
    'reported )'], round(pktnorx_wide_preamb(EbNo_idx)*100));
legend('wide beating (without preamble)',msg)
ylim([0 inf])

% Beating parameters
A2 = 1;                % Amplitude of transmitter 2, A1 is assumed to be 1
fbeat = 2/pduration;   % Beating frequency (Hz)

% Running the simulation
disp(['Calculating PER for 2 Concurrent Transmitters with narrow ' ...
    'beating and introducing the effect of the preamble...'])
[per_narrow_preamb, errindx_narrow_preamb,pktnorx_narrow_preamb] = ...
    ct_sim_run(M,Fs,nsamp,freqsep,A2,fbeat,EbNo,N,pktlen, preamblen);

% Plotting the results
figure
EbNo_idx = ceil(length(EbNo/2)); %Plot for an arbitrary EbNo index
plot(1:length(errindx_narrow(EbNo_idx,:)),errindx_narrow(EbNo_idx,:),'-*')
hold on;
plot(1:length(errindx_narrow_preamb(EbNo_idx,:)), ...
    errindx_narrow_preamb(EbNo_idx,:),'-*')
xlabel('Bit Index')
ylabel('Error Count')
msg = sprintf(['narrow beating (with preamble, %u%% of packets not ' ...
    'reported )'], round(pktnorx_narrow_preamb(EbNo_idx)*100));
legend('narrow beating (without preamble)', msg)
ylim([0 inf])

% As shown, there will be an important percentage of packets (though we are 
% simulating a worst-case scenario in which A2=1, energy diversity in real
% deployments would decrease the error rate) in which the preamble is 
% broken, and they are not even reported in the receiver as a packet with a 
% CRC error. They will be "phantom" transmissions. Therefore, even if some 
% level of error correction is applied at the application level, a robust 
% protocol based on CT should implement repetitions to assure a packet 
% reception is triggered and the payload is reported to the application.

%% Analyzing the effect of the bitrate
% The largest consequence of duplicating the bitrate is that the beating is
% widened, e.g. the packet experiences half the amount of valleys since the
% transmission takes half the time.

pktlen = 128;    % Packet length, in bits (with preamble)
preamblen = 32;  % Preamble length, in bits

% We arbitrarily fix the beating frequency to 20kHz
% Beating parameters
A2 = 1;                % Amplitude of transmitter 2, A1 is assumed to be 1
fbeat = 20e3;          % Beating frequency (Hz)

% Modulation parameters
M = 2;         % Modulation order
Rs = 1e6;      % Symbol rate (1Mbps)
freqsep = Rs;  % Frequency separation (optimal for noncoherent detection)
nsamp = 8;     % Number of samples per symbol
Fs = Rs*nsamp; % Sample rate (Hz)

% Running the simulation
disp(['Calculating PER for 2 Concurrent Transmitters at 1Mbps ' ...
    'and 20kHz beating...'])
[per_1M_20k, errindx_1M_20k,pktnorx_1M_20k] = ...
    ct_sim_run(M,Fs,nsamp,freqsep,A2,fbeat,EbNo,N,pktlen, preamblen);

% Modulation parameters
M = 2;         % Modulation order
Rs = 2e6;      % Symbol rate (2Mbps)
freqsep = Rs;  % Frequency separation (optimal for noncoherent detection)
nsamp = 8;     % Number of samples per symbol
Fs = Rs*nsamp; % Sample rate (Hz)

% Running the simulation
disp(['Calculating PER for 2 Concurrent Transmitters at 2Mbps ' ...
    'and 20kHz beating...'])
[per_2M_20k, errindx_2M_20k,pktnorx_2M_20k] = ...
    ct_sim_run(M,Fs,nsamp,freqsep,A2,fbeat,EbNo,N,pktlen, preamblen);

% Plotting the results
figure
EbNo_idx = ceil(length(EbNo/2)); %Plot for an arbitrary EbNo index
plot(1:length(errindx_1M_20k(EbNo_idx,:)),errindx_1M_20k(EbNo_idx,:),'-*')
hold on;
plot(1:length(errindx_2M_20k(EbNo_idx,:)), ...
    errindx_2M_20k(EbNo_idx,:),'-*')
xlabel('Bit Index')
ylabel('Error Count')
msg = sprintf(['20kHz beating at 1Mbps (%u%% of packets not ' ...
    'reported )'], round(pktnorx_1M_20k(EbNo_idx)*100));
msg2 = sprintf(['20kHz beating at 2Mbps (%u%% of packets not ' ...
    'reported )'], round(pktnorx_2M_20k(EbNo_idx)*100));
legend(msg, msg2)
ylim([0 inf])

% Which bitrate would perform better would actually depend on the frequency
% of the beating, which changes for every node pair (or it has a complex 
% shape for more than 2 CT). Differences in the sensitivity of the receiver
% at different bitrates should also be considered, i.e. a faster bitrate
% might have a shorter range due to the wider (noisier) spectrum.

%% Analyzing the effect of coding
% We will also modify the code of the simulator to analyze the effect of
% introducing convolutional coding with S=2 and S=8, inspired by the 500k
% and 125k modes of BLE5.
% Modulation parameters
M = 2;         % Modulation order
Rs = 1e6;      % Symbol rate (1Mbps)
freqsep = Rs;  % Frequency separation (optimal for noncoherent detection)
nsamp = 8;     % Number of samples per symbol
Fs = Rs*nsamp; % Sample rate (Hz)

% We arbitrarily fix the beating parmeters
A2 = 0.5;             % Amplitude of transmitter 2, A1 is assumed to be 1
fbeat = 2e3;          % Beating frequency (Hz)

% Running the simulation
disp(['Calculating PER for 2 Concurrent Transmitters at 500kbps ' ...
    'and 2kHz beating with coding S=2...'])
[per_500k_2k, errindx_500k_2k,pktnorx_500k_2k] = ...
    ct_sim_run_fec(M,Fs,nsamp,freqsep,A2,fbeat,EbNo,N,pktlen, ...
    preamblen,false);

% Running the simulation
disp(['Calculating PER for 2 Concurrent Transmitters at 125kbps ' ...
    'and 2kHz beating with coding S=8...'])
[per_125k_2k, errindx_125k_2k,pktnorx_125k_2k] = ...
    ct_sim_run_fec(M,Fs,nsamp,freqsep,A2,fbeat,EbNo,N,pktlen, ...
    preamblen,true);

% Plotting the results
figure
plot(EbNo,per_500k_2k,'-*')
hold on
plot(EbNo,per_125k_2k,'-*')
set(gca,'YScale','log')
xlabel('Eb/N0')
ylabel('PER')
legend('PER 500kbps, S=2 (2kHz beating, A2=0.5)', ...
    'PER 125kbps, S=8 (2kHz beating, A2=0.5)')

%% Simulator code
function [per, errindx, pktnorx] = ct_sim_run(M,Fs,nsamp,freqsep,A2, ...
    fbeat,EbNo,N,pktlen,preamblen)
    arguments
        M (1,1)
        Fs (1,1)
        nsamp (1,1)
        freqsep (1,1)
        A2 (1,1)
        fbeat (1,1)
        EbNo (1,:)
        N (1,1)
        pktlen (1,1)
        preamblen (1,1) {mustBeLessThanOrEqual(preamblen,pktlen)}
    end

    pktnorx = zeros(1,length(EbNo));
    errindx = zeros(length(EbNo),pktlen);
    t = 0:1/Fs:(N*nsamp/Fs-1/Fs);
    per = zeros(1,length(EbNo));
    data = randi([0 M-1],N,1);
    sig_beating = fskmod(data,M,freqsep,nsamp,Fs);

    for i = 1:length(sig_beating)
        % We consider a different (and random) initial beating phase for 
        % each packet transmission
        if mod(i-1,pktlen*nsamp) == 0
            phase = 2*pi*rand(); 
        end
        % We introduce the beating distortion to the modulated signal
        sig_beating(i) = sig_beating(i)*...
            (1-A2+2*A2*abs(cos(2*pi*((fbeat/2)*t(i))+phase)));
    end

    for i = 1:length(EbNo)
        % We add AWGN noise to the concurrent transmission
        sig_beating_noisy  = awgn(sig_beating,EbNo(i)+10*log10(log2(M))...
            -10*log10(nsamp),0,'dB');

        dataOut = fskdemod(sig_beating_noisy,M,freqsep,nsamp,Fs);

        total_pkt = 0;
        for j = 1:pktlen:(N-pktlen+1)
            first_pkt_error = true;
            total_pkt = total_pkt+1;
            for s = j:(j+pktlen-1)
                if data(s) ~= dataOut(s)
                    if s-j+1 < preamblen
                        per(i) = per(i)+1;
                        pktnorx(i) = pktnorx(i)+1;
                        break;
                    else
                        errindx(i,s-j+1) = errindx(i,s-j+1)+1; 
                        if first_pkt_error
                            per(i) = per(i)+1;
                            first_pkt_error = false;
                        end
                    end
                end
            end
        end
        pktnorx(i) = pktnorx(i)./total_pkt;
        per(i) = per(i)./(N/pktlen);
    end
end
%% Variation for coding
function [per, errindx, pktnorx] = ct_sim_run_fec(M,Fs,nsamp,freqsep, ...
    A2,fbeat,EbNo,N,pktlen,preamblen,S8)
    arguments
        M (1,1)
        Fs (1,1)
        nsamp (1,1)
        freqsep (1,1)
        A2 (1,1)
        fbeat (1,1)
        EbNo (1,:)
        N (1,1)
        pktlen (1,1)
        preamblen (1,1) {mustBeLessThanOrEqual(preamblen,pktlen)}
        S8 (1,1) logical
    end

    if S8
        coderate = 1/8;
    else
        coderate = 1/2;
    end

    pktnorx = zeros(1,length(EbNo));
    errindx = zeros(length(EbNo),pktlen);
    t = 0:1/Fs:(N*nsamp*(1/coderate)/Fs-1/Fs);
    per = zeros(1,length(EbNo));
    data = randi([0 M-1],N,1);
    tPoly = poly2trellis(4,[17 15]);
    tx_data = zeros(1,(1/coderate)*length(data));
    coded_data = zeros(1,2*length(data));

    % We code the bitstream with a 1:2 convolutional encoder (S=2)
    for j = 1:pktlen:(N-pktlen+1)
        coded_data(2*(j-1)+1:2*(j+pktlen-1)) = convenc( ...
            data(j:j+pktlen-1),tPoly);
    end
    
    % If we choose S=8 we add a simple orthogonal mapping of 1:4 ratio 
    if S8
        for i=1:length(coded_data)
            if coded_data(i)==0
                tx_data((4*(i-1)+1):(4*(i-1)+1+3))= [0 0 1 1];
            else
                tx_data((4*(i-1)+1):(4*(i-1)+1+3))= [1 1 0 0];
            end
        end         
    else
        tx_data = coded_data;
    end

    sig_beating = fskmod(tx_data,M,freqsep,nsamp,Fs);

    for i = 1:length(sig_beating)
        % We consider a different (and random) initial beating phase for 
        % each packet transmission
        if mod(i-1,pktlen*nsamp*(1/coderate)) == 0
            phase = 2*pi*rand(); 
        end
        % We introduce the beating distortion to the modulated signal
        sig_beating(i) = sig_beating(i)*...
            (1-A2+2*A2*abs(cos(2*pi*((fbeat/2)*t(i))+phase)));
    end

    for i = 1:length(EbNo)
        % We add AWGN noise to the concurrent transmission
        sig_beating_noisy  = awgn(sig_beating,EbNo(i)+10*log10(log2(M))...
            -10*log10(nsamp),0,'dB');

        dataOut = fskdemod(sig_beating_noisy,M,freqsep,nsamp,Fs);

        % If S=8 we perform the 4:1 demapping 
        if S8
            for z=1:length(coded_data)
                if sum(dataOut((4*(z-1)+1):(4*(z-1)+1+3)).* ...
                        [0 0 1 1]) > sum(dataOut((4*(z-1)+1): ...
                        (4*(z-1)+1+3)).*[1 1 0 0])
                    dataOut(z) = 0;
                elseif sum(dataOut((4*(z-1)+1):(4*(z-1)+1+3)).* ...
                        [0 0 1 1]) < sum(dataOut((4*(z-1)+1): ...
                        (4*(z-1)+1+3)).*[1 1 0 0])
                    dataOut(z) = 1;
                else
                    dataOut(z) = randi([0, 1], [1, 1]);
                end 
            end
        end

        total_pkt = 0;
        for j = 1:pktlen:(N-pktlen+1)

            % We decode the received bitstream with a Viterbi decoder 2:1
            dataOut_dec = vitdec(dataOut(2*(j-1)+1:2*(j+pktlen-1)), ...
                tPoly,5*(4-1),'trunc','hard');
            
            first_pkt_error = true;
            total_pkt = total_pkt+1;
            for s = j:(j+pktlen-1)
                if data(s) ~= dataOut_dec(s-j+1)
                    if s-j+1 < preamblen
                        per(i) = per(i)+1;
                        pktnorx(i) = pktnorx(i)+1;
                        break;
                    else
                        errindx(i,s-j+1) = errindx(i,s-j+1)+1; 
                        if first_pkt_error
                            per(i) = per(i)+1;
                            first_pkt_error = false;
                        end
                    end
                end
            end
        end
        pktnorx(i) = pktnorx(i)./total_pkt;
        per(i) = per(i)./(N/pktlen);
    end
end