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
% amplitude modulation:
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

%%%
% First of all, we perform a basic sanity check, by setting to zero the
% amplitude of the second transmitter, effectively removing it from the
% simulation, and fixing the packet length to 1 bit. In that case, the
% result should be the well-known Bit Error Rate (BER) for non-coherent FSK
% receivers:
%
% $$ BER_{BFSK} =\frac{1}{2}\exp\left(-\frac{E_{b}}{2N_{0}}\right) $$

% Beating parameters
A2 = 0;        % Amplitude of transmitter 2, A1 is assumed to be 1
fbeat = 0.01;  % Beating frequency (Hz)

% Simulation parameters
EbNo = 0:1:12; % Eb/No (dB), desired simulation range 
N = 200000;   % Transmitted Bits
p_length = 1;  % Packet length, in bits

% Running the simulation
per = ct_sim_run(M,Fs,nsamp,freqsep,A2,fbeat,EbNo,N,p_length);

% Analytical results for comparison
EbNo_linear = 10.^(EbNo/10);
ber_theory = 0.5*exp(-EbNo_linear/2);

% Plotting the results
plot(EbNo,per,'b*')
hold on;
plot(EbNo,ber_theory,'b-')
set(gca, 'YScale', 'log')
xlabel('Eb/N0')
ylabel('BER')
legend('Simulated BFSK BER', 'Analytical BFSK BER')

%% Simulator code
function per = ct_sim_run(M,Fs,nsamp,freqsep,A2,fbeat,EbNo,N,...
    p_length)
    t=0:1/nsamp:(N-1/nsamp);
    per=zeros(1,length(EbNo));
    data = randi([0 M-1],N,1);
    sig_beating = fskmod(data,M,freqsep,nsamp,Fs);

    for i=1:length(sig_beating)
        % We consider a different (and random) intiial beating phase for 
        % each packet transmission
        if mod(i-1,p_length*nsamp) == 0
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

        for j=1:p_length:(N-p_length)
            for s=j:(j+p_length-1)
                if data(s) ~= dataOut(s)
                    per(i) = per(i)+1;
                    break;
                end
            end
        end
        per(i) = per(i)./(N/p_length);
    end
end  