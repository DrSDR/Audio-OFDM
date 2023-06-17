

% ofdm decoder
clear all
close all
clc
plotcrap = 1;

FILTERSPEC = '.wav';
TITLE = 'Pick a Recorded OFDM IQ wav file';
FILE = 'D:\gnu_radio_work';

[FILENAME, PATHNAME, FILTERINDEX] = uigetfile(FILTERSPEC, TITLE, FILE);
xstring = [PATHNAME FILENAME];

[x,fswav] = audioread(xstring);

% tx bandwidth

fs = 8e3;  % sample rate of ofdm signal 
% pulse at start for sync 
N = 8000;  %number of samples for chirp sync pulse
t = [0:7999]/fs;
f0 = 300;
f1 = 3000;
T = 1;
c = (f1 - f0)/T;
x1 = (c/2)*t.^2 + f0*t;
pream = sin(2*pi*x1);

Nguard = 512;   % time between pulse and preamble fft data

Nchar = 256;  % make power of 2
bits = Nchar * 8;
N2 = bits/2;
Nfft = bits + 1024;
Nfft = (2*Nfft)  - 1 ; 




sigtime = Nguard + 2*Nfft + 2*Nfft + 2*Nguard;





x = x(:,1);
x = x.';






if fswav ~= fs
    n = gcd(fs,fswav);
    p = fs / n;
    q = fswav / n;
    x = resample(x,p,q);   % resample file to expected sample rate
end


hpw = conj(pream(end:-1:1));
xdet = filter(hpw,1,x);

if plotcrap
    figure(100)
    plot(xdet)
    title('chirp detection')
end

[maxv, maxi] = max(abs(xdet));
xstart = maxi ;
xend = xstart + sigtime ;

if plotcrap
    figure(101)
    xplot = abs(x) / max(abs(x));
    xplot = 20 * log10(xplot);
    yplot = length(xdet);
    yplot = zeros(1,yplot);
    yplot(maxi) = -40;
    plot(xplot)
    hold on
    plot(yplot)
    hold off
end





x = x(xstart:xend);
a = Nguard  + round(0.8*Nfft) + 1 ;
b = a + Nfft   - 1 ;
xp = x(a:b);

a1 =  Nguard + 2*Nfft + round(0.8*Nfft) + 1 ;
b2 = a1 + Nfft  - 1 ;
xd = x(a1:b2);







Xp = fft(xp);
Xp = Xp( 1:floor(end/2) );
Xd = fft(xd);
Xd = Xd( 1:floor(end/2) );

if plotcrap
    n = length(Xd);
    fvec = linspace(0,fs/2,n);
    figure(103)
    plot(fvec,20*log10(abs(Xp)))
    hold on
    plot(fvec,20*log10(abs(Xd)))
    legend('Pilot FFT', 'Data FFT');
    hold off
end










det = angle(Xp ./ Xd);          % ratio to get delta phase diff
det = abs(det);

thres = pi/2;

% det = ( det <= thres) .* 1   +  (det > thres) .* 0;
det = det(513:end);
det = det(1:bits);   % bits for char  

if plotcrap
    figure(104)
    plot(det)
    title('demodulate phase vector')
end





det(det <= thres) = 1;
det(det > thres) = 0;

if plotcrap
    figure(105)
    plot(det)
    title('decoded bits vector')
end


det = reshape(det,8,[]);
det = det';


charmes = zeros(1,Nchar);
w = [7:-1:0];
w = 2.^w;


for k = 1:Nchar
    x = det(k,:);
    x = x .* w;
    x = sum(x);
    charmes(k) = x;
end



    clc

   xstr =  char(charmes)

   msgbox(xstr,'replace')
















