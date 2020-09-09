%This file is for experimenting with Fourier Series and Fourier Transform
%implemtnations. This code turned into the Fourier Series C++ code. 

Fs =  100;           % Sampling frequency                    
T = 1/Fs;             % Sampling period      
L = 200;             % Length of signal
t = (0:L-1)*T;        % Time vector#
d = Fs/L;
f = (0:L-1)*d; % <-- Don't understand this 
u = 2*pi*f;
%noise = 1.2*randn(1,L);
y = 4*sin(2*pi*4*t+1)+sin(2*pi*3*t+0.3)+cos(2*pi*9*t+0.2)*2;

%Amps:5,1,7. Freq: 4,6,3
yf = fft(y)/L;

[out,idx] = sort(yf,'descend');
yf_sorted=out;
%u=u(idx);

yy = complex(zeros(size(y)));
terms = 4;
%Complex Fourier series with matrix operations.
% for k = 1:terms
%     yy = yy + yf_sorted(k)*exp(1i*u(k)*t);
% end

% Complex Fourier series with no matrix operations.
for k = 1:terms
    pr = 2 * pi * ((idx(k)-1) * (Fs / L));
    for i = 1:L
        t_ = 1.0 / Fs * (i-1);
        debug2 = 1i*pr*t_;
        yy(i) = yy(i) + yf_sorted(k)*exp(debug2);
    end
end

% for i = 1:terms
%     k=idx(i);
%     pr = 2 * pi * (k * (Fs / L));
%     for x = 1:L
%         t = 1.0 / Fs * x;
%         yy(i) = yf_sorted(i) * exp(1i * pr * t);
%     end
% end

%plot(t,yy,t,y);

%plot(t,y,t,yy)     % add offset for visual purposes

f=figure('visible','on');

tiledlayout(3,1);
nexttile
plot(t,y);
title("original function");
nexttile
soooslow=2*abs(yf(1:L/2+1));
plot(Fs*(0:(L/2))/L,soooslow);
title("Frequency spectrum");
nexttile
plot(t,real(yy));
title("Synthesis from spectrum ("+terms+ " terms)");
% error between synthesis and original.
disp(corrcoef(y,real(yy)));
% create datset from y and t
vals = y.';
time = t.';
data=horzcat(time,vals);

