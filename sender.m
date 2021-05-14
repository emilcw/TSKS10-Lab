function x = sender(xI, xQ)
% SENDER     Upsample xI, xQ, apply chirp and generate 
%            signal ready to be sent to channel.
%
%   SENDER(xI, xQ) returns x(t) upsampled with an applied 
%                  pulseform and the beginning
%                  so that it can be sent over a channel 
%                  via I/Q- modulation.
%    
%-------------------------Internal variables-------------------------------
f1 = 85e3;
f2 = 105e3;
bandwidth = 5e3;
fs = 400e3;
fn = fs / 2;
Ts = 1/fs;

%------------------------ Upsample to 400kHz-------------------------------
sample_factor = 20;
% Multiply by sample_factor to compensate for loss of signal energy
xI_upsamp = upsample(xI, sample_factor) * sample_factor;
xQ_upsamp = upsample(xQ, sample_factor) * sample_factor;

% Lowpass-filter to make it interpolate
n_lp = 100; %According to task 
W = (2 * bandwidth) / (fn);
[b,a] = fir1(n_lp,W,'low');
xI_low = filter(b,a,xI_upsamp);
xQ_low = filter(b,a,xQ_upsamp);

%Compensate for filter delay
%xI_low = xI_low((n_lp/2)+1:end);
%xQ_low = xQ_low((n_lp/2)+1:end);

%----------------------Apply pulseform ------------------------------------

%Generate chirp-signal
t1 = 0:Ts:1-Ts;
f0 = 1;         %Start frequency
epsilon = 10;    %Chirp-rate
chirp = cos(2*pi*f0*(1+(epsilon.*t1)).*t1);
chirp_transposed = transpose(chirp);

% Apply chirp to xI_low and zeroes to xQ_low to match dimensions.
xI_with_chirp = [chirp_transposed; xI_low];
xQ_with_zeros = [zeros(length(chirp_transposed),1); xQ_low];

%------------------Generate signal to send via I/Q-------------------------

%I/Q - modulation
t2 = 0:Ts:6-Ts;
fc = (f1 + f2) / 2;
I_carry = transpose(cos(2*pi*fc*t2));
Q_carry = transpose(-sin(2*pi*fc*t2));
xI_send = xI_with_chirp .* I_carry;
xQ_send = xQ_with_zeros .* Q_carry;

x = xI_send + xQ_send;