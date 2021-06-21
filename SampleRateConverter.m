% Sample Rate Converter using farrow filter structure
%
% References:
% [1] https://ccrma.stanford.edu/~jos/pasp/Farrow_Structure.html
% [2] https://ccrma.stanford.edu/~jos/pasp/Lagrange_Interpolation.html
% [3] http://users.spa.aalto.fi/vpv/publications/vesa_phd.html
% [4] https://www.mathworks.com/help/dsp/ref/dsp.farrowrateconverter-system-object.html

clc; clear all; close all;


%% Input signal

fc = 1.5e6;
fs_ip = 20e6; % Input sample rate
fs_op = 21e6; % Output(desired) sample rate

[L, M] = rat(fs_op/fs_ip);

t_str = 0;
t_end = 10e-6; 

t_ip = t_str:1/fs_ip:(t_end - 1/fs_ip);
t_op = t_str:1/fs_op:(t_end - 1/fs_op);

signal = cos(2*pi*fc*t_ip);

figure('units','normalized','outerposition',[0 0 1 1]);
plot(t_ip,signal,'b',t_ip,signal,'ro');
delta = 0.1;
hold on;
stem(t_ip,signal,'r');
hold off;
ylim([min(signal)-delta max(signal)+delta]);
title('Input signal'); legend('Continuous','Sampled (20MHz)'); xlabel('Time'); ylabel('Amplitude');


%% Farrow filter Design

order = 3;
coeffs = designFarrowFilt(order); % Native Implementation

farrowRateConv_3rd = dsp.FarrowRateConverter('InputSampleRate',fs_ip, ...
    'OutputSampleRate',fs_op,'PolynomialOrder',3); % MATLAB reference design
[L_ref, M_ref, C_ref, filtObj] = farrowRateConv_3rd.designPolynomialSRC(fs_ip, fs_op,0,1,order); 

tol = 1e-8;
if (sum(sum(coeffs - C_ref,1),2) < tol)
    disp('Filter design has been validated!');
end

%% Applying Farrow filter to Input signal to get rate converted signal

signalRC = applyFilter(signal, coeffs, L, M); % Native Implementation

signalRC_ref = farrowRateConv_3rd(signal.'); % MATLAB reference 

if (sum(sum(signalRC - signalRC_ref,1),2) < tol)
    disp('Filter output has been validated!');
end

%% Frequency Response of farrow filter

W = linspace(0,fs_op*3,2048);  % Define the frequency range analysis
Fs1 = fs_ip*33;  % The equivalent single stage filter is clocked at 3.53 MHz
hfvt = fvtool(farrowRateConv_3rd,'FrequencyRange','Specify freq. vector', ...
    'FrequencyVector',W,'Fs',[Fs1],'NormalizeMagnitudeto1','on','Color','white');
legend(hfvt, '3rd-Order Farrow Interpolator', 'Location','NorthEast');

%% Group delay of farrow filter

grpdelay(farrowRateConv_3rd,128,fs_op);
[D,F] = grpdelay(farrowRateConv_3rd,128,21e6);
delay = round(D(end));
delay = delay+2; % adjustment to syncronize with signal

%% Group delay compensation

t_op = t_op(1:end-delay);
signalRC(1:delay) = [];
signalRC_ref(1:delay) = [];

%% Resampled signal(after group delay compensation) plotted against Original signal

figure('units','normalized','outerposition',[0 0 1 1]);
plot(t_ip,signal,'b');
delta = 0.1;
hold on;
stem(t_ip,signal,'r');
hold on;
stem(t_op,signalRC,'g');
hold off;
ylim([min(signal)-delta max(signal)+delta]);
title('Original signal vs resampled signal');
legend('Signal(continuous)','Sampled signal (20MHz)','Farrow rate converted (21MHz)');
title('Input signal and Resampled signal after compensating for group delay');
xlabel('Time(Sec)'); ylabel('Amplitude');

%% Farrow Filter Design
% References:
% [5] http://users.spa.aalto.fi/vpv/publications/vesan_vaitos/ch3_pt2_lagrange.pdf

function coeffs = designFarrowFilt(order)

N = order;

U = zeros(N+1,N+1); % Eq (3.94)
for i=1:size(U,1)
    for j=1:size(U,2)
        U(i,j) = (i-1)^(j-1);
    end
end

Q = pinv(U); % Eq (3.113)

T = zeros(N+1,N+1); % Eq (3.109)
for i=1:size(T,1)
    for j=1:size(T,2)
        
        n = j-1;
        m = i-1;       
        if n >= m
           T(i,j) = floor(N/2)^(n-m)*nchoosek(n, m);
        else
            T(i,j) = 0;
        end
            
    end
end

Q = T*Q;
coeffs = flipud(Q);

end

%% Function that applies farrow filter to signal

function signalOut = applyFilter(signalIn, coeffs, L, M)

C = coeffs;
R = zeros(size(coeffs,1),1);
sigLen = length(signalIn);

if size(signalIn,2)>size(signalIn,1)
    signalIn = signalIn.';
end

InterpStep = L;
Linv   = 1/L;
Np = size(coeffs,1)-1;
outIdx = 1;

for ipIdx = 1:sigLen
    
    % Update internal states of filter with incoming signal
    R = [signalIn(ipIdx); R(1:end-1)];
    
    % Apply farrow filter coefficients and update current state of filter
    RC = C * R;
    
    % Fractional Delay loop processing
    while (InterpStep > 0)
        % Compute fractional delay step
        fd = InterpStep .* Linv;
        
        % Multiply state by fractional delay step
        signalOut(outIdx) = fd * RC(1);
        for pIdx = 2:Np
            signalOut(outIdx) = fd * (RC(pIdx) + signalOut(outIdx));
        end
        signalOut(outIdx) = signalOut(outIdx) + RC(Np+1);

        % Update Interpolation step and continue to next step
        InterpStep = InterpStep - M;
        outIdx  = outIdx + 1;
    end
      
    InterpStep = InterpStep + L;
end

end
