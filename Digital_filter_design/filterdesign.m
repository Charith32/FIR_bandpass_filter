%filter specifications
close all;
clear all;
clc;
Ap = 0.08;  %maximum passband ripple
Aa = 43;    %minimum stopband attenuation
wp1 = 400;  %lower passband edge
wp2 = 700;  %upper passband edge
wa1 = 250;  %lower stopband edge
wa2 = 800;  %upper stopband edge
ws = 2400;  %sampling frequency

Bt = min((wp1-wa1),(wa2-wp2));
wc1 = wp1 - Bt/2;
wc2 = wp2 + Bt/2;
T = 2*pi/ws; %sampling period

deltaph = (10^(0.05*Ap) - 1)/(10^(0.05*Ap) + 1);
deltaah = 10^(-0.05*Aa);
delta = min(deltaph,deltaah);

Aaact = -20*log10(delta); %actual stopband attenuation

%choosing alpha parameter
if Aaact <= 21
    alpha = 0;
elseif Aaact <= 50
    alpha = 0.5842*((Aaact - 21)^0.4) + 0.07886* (Aaact - 21);
else
    alpha = 0.1102*(Aaact - 8.7);
end

%choosing D parameter
if Aaact <= 21
    D = 0.9222;
elseif Aaact > 21
    D = (Aaact - 7.95)/14.36;
end    
% Choosing N value
val = ((ws*D)/Bt) + 1 ;
val = ceil(val);
if mod(val,2) == 0 
    N = val + 1;
else
    N = val;
end

%filter n values
n = -(N-1)/2 : 1 : (N-1)/2;
beta = alpha * sqrt(1 - ((2*n)/(N-1)).^2);
%Io(beta) and Io(alpha)
Ioalpha = bessfun(alpha);
Iobeta = bessfun(beta);
%generating window function
wknt = Iobeta/Ioalpha;

%Impulse response representation of kaizer window
figure;
stem(n,wknt);
xlabel('n')
ylabel('Wk(nT) - Amplitude')
title('Kiaser window - Impulse response')

%generating the ideal frequency response
hnt0 = (2/ws)*(wc2-wc1);
rrange = 1 : (N-1)/2;
hntright = (1./(rrange*pi)).*(sin(wc2.*rrange*T) - sin(wc1.*rrange*T)) ;
lrange = -(N-1)/2 : -1;
hntleft = (1./(lrange*pi)).*(sin(wc2.*lrange*T) - sin(wc1.*lrange*T)) ;
hnt = [hntleft,hnt0,hntright];
figure;
stem(n,hnt);
xlabel('n');
ylabel('hi[n]');
title('Time domain of the ideal band pass filter');

% combining the ideal bandpass filter with kaiser window to get the h[n]
hn = hnt.*wknt;
figure;
n = 1:N;
stem(n,hn)
xlabel('n')
ylabel('h[n]')
title('Time domain of the filter transfer function(h[n]) using windowing function')

%obtaining the magnitude response of h[n] 
[hnMag,f] = freqz(hn);
w = f/T;
loghnMag = 20*log10(abs(hnMag));
figure;
plot(w,loghnMag);
grid on;
xlabel('frequency (rad/s)');ylabel('Magnitude(dB)')
title('Freq domain of the filter')

%obtaining the magnitude response for frequncies in the passband
figure;
plot(w,loghnMag);
grid on;
axis([wc1 wc2 -0.1 0.1]);
xlabel('frequency (rad/s)'); ylabel('Magnitude(dB)');
title('Magnitude response in the passband of filter')

%Generating the input signal
w1 = wc1/2;
w2 = wc1 + (wc2-wc1)/2;
w3 = wc2 + (ws/2 - wc2)/2;
L = 2048;   %sample size
n = 0:L-1;
f = (-L/2:L/2-1)*ws/L;
x_n=sin(n*w1*T)+sin(n*w2*T)+sin(n*w3*T);   %generated input
t_n=sin(n*w2*T);    %ideal output signal


X_N = fftshift(fft(x_n,L));
T_N = fftshift(fft(t_n,L));
H_N = fftshift(fft(hn,L));

Y_N = X_N.*H_N; % frequency domain multiplication of the signals
%convolution of x_n and hn
m=length(x_n);
n=length(hn);
X=[x_n,zeros(1,n)];
H=[hn,zeros(1,m)];
Y = zeros(n+m-2);
for i=1:n+m-1
    Y(i)=0;
    for j=1:m
        if(i-j+1>0)
            Y(i)=Y(i)+X(j).*H(i-j+1);
        else
        end
    end
end
y_n = Y;


figure;
subplot(3,1,1);plot(n,x_n);
ylim([-2 2]);xlim([0 L]);
xlabel('Samples (n)');ylabel('Amplitude');
title('Time domain of the input signal');
subplot(3,1,2);plot(n,y_n((N-1)/2:L+(N-1)/2-1));
ylim([-2 2]);xlim([0 L]);
xlabel('Samples (n)');ylabel('Amplitude');
title('Time domain output of the FIR BPF');
subplot(3,1,3);plot(n,t_n);
ylim([-2 2]);xlim([0 L]);
xlabel('Samples (n)');ylabel('Amplitude');
title('Time domain of the ideal signal(Expected)');


figure;
subplot(3,1,1);plot(f,abs(X_N/L));
xlabel('Frequency (rad/s)');
ylabel('Magnitude');
title('Frequency domain of the input signal');
subplot(3,1,2);plot(f,abs(Y_N/L));
xlabel('Frequency (rad/s)');
ylabel('Magnitude');
title('Frequency domain output of the FIR BPF');
subplot(3,1,3);plot(f,abs(T_N/L));
xlabel('Frequency (rad/s)');
ylabel('Magnitude');
title('Frequency domain of the ideal signal(Expected)');







































    
    











    
    
    













