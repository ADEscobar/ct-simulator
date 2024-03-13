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
% consider non-coherent reception of binary FSK transmissions. 

% Modulation parameters
M = 2;         % Modulation order
Fs = 40;       % Sample rate (Hz)
nsamp = 8;     % Number of samples per symbol
freqsep = 5;   % Frequency separation (Hz)

% Simulation parameters
EbNo = 0:2:12; % Eb/No (dB), desired simulation range 
N = 100000;    % Transmitted Bits
plength = 1;   % Packet length, in bits

%%%
% First of all, we perform a basic sanity check, by setting to zero the
% amplitude of the second transmitter, effectively removing it from the
% simulation, and fixing the packet length to 1 bit. In that case, the
% result should be the well-known Bit Error Rate (BER) for non-coherent FSK
% receivers:
%
% $$ BER_{BFSK}=\frac{1}{2}\exp\left(-\frac{E_{b}}{2N_{0}}\right) $$

% Beating parameters
A2 = 0;        % Amplitude of transmitter 2, A1 is assumed to be 1
fbeat = 0.01;  % Beating frequency (Hz)

% Running the simulation
disp('Calculating BER for 1 Transmitter...')
ber_1 = ct_sim_run(M,Fs,nsamp,freqsep,A2,fbeat,EbNo,N,plength);

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
fbeat = 0.01;  % Beating frequency (Hz)

% Running the simulation
disp('Calculating BER for 2 Concurrent Transmitters...')
ber_2 = ct_sim_run(M,Fs,nsamp,freqsep,A2,fbeat,EbNo,N,plength);

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
fbeat = 0.01;  % Beating frequency (Hz)

% Running the simulation
disp('Calculating BER for 2 Concurrent Transmitters with energy delta...')
ber_2_en_delta = ct_sim_run(M,Fs,nsamp,freqsep,A2,fbeat,EbNo,N,plength);

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
% beating period four times shorter for narrow beating. 

plength = 128;                % Packet length, in bits
pduration = plength*nsamp/Fs; % Packet duration in the air, in seconds

% Beating parameters
A2 = 1;                % Amplitude of transmitter 2, A1 is assumed to be 1
fbeat = 0.5/pduration; % Beating frequency (Hz)

% Running the simulation
disp('Calculating PER for 2 Concurrent Transmitters with wide beating...')
per_2_wide = ct_sim_run(M,Fs,nsamp,freqsep,A2,fbeat,EbNo,N,plength);

% Beating parameters
A2 = 1;                % Amplitude of transmitter 2, A1 is assumed to be 1
fbeat = 4/pduration;   % Beating frequency (Hz)

% Running the simulation
disp(['Calculating PER for 2 Concurrent Transmitters ' ...
    'with narrow beating...'])
per_2_narrow = ct_sim_run(M,Fs,nsamp,freqsep,A2,fbeat,EbNo,N,plength);

% Plotting the results
figure
plot(EbNo,per_2_wide,'-*')
hold on;
plot(EbNo,per_2_narrow,'-*')
set(gca,'YScale','log')
xlabel('Eb/N0')
ylabel('BER')
legend('BFSK PER (wide beating)','BFSK PER (narrow beating)')

%% Simulator code
function per = ct_sim_run(M,Fs,nsamp,freqsep,A2,fbeat,EbNo,N,plength)
    t=0:1/nsamp:(N-1/nsamp);
    per=zeros(1,length(EbNo));
    data = randi([0 M-1],N,1);
    sig_beating = fskmod(data,M,freqsep,nsamp,Fs);

    for i=1:length(sig_beating)
        % We consider a different (and random) initial beating phase for 
        % each packet transmission
        if mod(i-1,plength*nsamp) == 0
            phase = 2*pi*rand(); 
        end
        % We introducethe beating distortion to the modulated signal
        sig_beating(i) = sig_beating(i)*...
            (1-A2+2*A2*abs(cos(2*pi*((fbeat/Fs)/2)*t(i)+phase)));
    end

    for i=1:length(EbNo)
        % We add AWGN noise to the concurrent transmission
        sig_beating_noisy  = awgn(sig_beating,EbNo(i)+10*log10(log2(M))...
            -10*log10(nsamp),0,'dB');

        dataOut = fskdemod(sig_beating_noisy,M,freqsep,nsamp,Fs);

        for j=1:plength:(N-plength)
            for s=j:(j+plength-1)
                if data(s) ~= dataOut(s)
                    per(i) = per(i)+1;
                    break;
                end
            end
        end
        per(i) = per(i)./(N/plength);
    end
end  