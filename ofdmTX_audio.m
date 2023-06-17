

clear all
close all

fs = 8e3;  % sample rate of ofdm signal 
% pulse at start for sync 
N = 8000;  %number of samples for chirp sync pulse
t = [0:7999]/fs;
f0 = 300;
f1 = 3000;
T = 1;
c = (f1 - f0)/T;
x = (c/2)*t.^2 + f0*t;
pream = sin(2*pi*x);
pream = pream / max(pream);






% samples from sync pulse to ofdm data
guardN = 512;


td = 2;
tdsamples = round(td * fs);
tdvec = zeros(1,tdsamples);


Nchar = 256;  % make power of 2,  number of char for text
bits = Nchar * 8;
N2 = bits/2;
Nfft = bits + 1024;  % bits plus guard freq bins

% get message from txt file
fid = fopen('C:\SDR_Work\textfile.txt');
mtxt = fread(fid);
fclose(fid);
mtxt = mtxt';
if length(mtxt) >= Nchar
    mtxt = mtxt(1:Nchar);
else
    z = length(mtxt);
    z = Nchar - z;
    mtxt = [mtxt zeros(1,z)];
end
mbits = dec2bin(mtxt,8);
mbits = reshape(mbits',1,[]);
mbits = mbits(1:bits);
txbits = zeros(1,bits);
for k = 1:bits
    if mbits(k) == '1'
        txbits(k) = 1;
    else
        txbits(k) = 0;
    end
end




% make pilot vector
xp = randn(1,Nfft);
xp = xp / max(xp);
xp = xp * 12;


% pilot vector with 512 guard on each side
xp = exp(1i*xp);
xp(1:512) = 0;
xp(end-511:end) = 0;

% make data vector bpsk,  1 = 1,  0 = -1
xd = txbits;  % put binary vector here N = 2048 bits
xd(xd == 0) = -1;
xd = [zeros(1,512)   xd    zeros(1,512) ];  %512 guard on each side






% relate data vector to pilot vector
xd = xp .* xd;  % e^a * e^b =  e^(a + b)




xp1 = conj( xp(end:-1:2));
xp = [xp1 xp];

xd1 = conj(xd(end:-1:2));
xd = [xd1 xd];


% pilot time vector
xp = ifftshift(xp);
xp = ifft(xp);
xp = xp / max(abs(xp));

% data time vector
xd = ifftshift(xd);
xd = ifft(xd);
xd = xd / max(abs(xd));

% put the whole movie together
xz = zeros(1,guardN);
xt = [tdvec pream xz xp xp xd xd tdvec];

% ensure bandwidth, reduce out of band signals
hlpf = fir1(64,0.9);
xt = filter(hlpf,1,xt);
xt = xt / max(abs(xt));  % 1xN

datafile = [ real(xt).' ];



 audiowrite('c:\SDR_Work\ofdmtest256Char8khzAudio.wav',datafile,fs,'BitsPerSample',16);
% 



figure(12)
plot(real(xt))
hold on
plot(imag(xt))
hold off
% 
% figure(13)
% plot(20*log10(abs(xt)))
% 
% 
% 
% 
Xp = fft(xp);
Xp = Xp( 1:floor(end/2) );

figure(56)
plot(abs(Xp))

Xd = fft(xd);
Xd = Xd( 1:floor(end/2) );


det = angle(Xp ./ Xd );

 det = abs(det);
n = length(det);
f = linspace(0,fs/2,n);

thres = pi/2;

% det = ( det >= thres) .* 1   +  (det < thres) .* 0;


figure(34)
plot(f,det)
title('phase results')


bits = det(513:end-511);

figure()
plot(bits)








