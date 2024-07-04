N=4096;
F=1000;
T=1/F;

% creating tx signal
bits = rand(N,1)>0.5;
bitT1 = bits(1:4:N);
bitT2 = bits(2:4:N);
bitT3 = bits(3:4:N);
bitT4 = bits(4:4:N);

Srs = 2*bitT1-1; %-2*bitT1+1
Sra = 2*bitT2+1; % 2*bitT2-3
Sis = 2*bitT3-1;
Sia = 2*bitT4+1;
Si = Srs.*Sra;
Sq = Sis.*Sia;
S = (Si+1i*Sq)/sqrt(10);


Ns=size(S,1);
Ncp=round(Ns/10);
Stx = ifft(fftshift(S));
Scp = Stx(Ns-Ncp+1:Ns);
NewS = cat(1,Scp,Stx);


% noise
snr = [70:-0.5:-20];
N0 = 10.^(-snr./10);
noise = sqrt((N0)/2).*randn(Ns+Ncp,1)+1i*sqrt((N0)/2).*randn(Ns+Ncp,1);

NewS = NewS+noise;

% processing rx signal
Srx = NewS(size(Scp,1)+1:size(Stx,1)+size(Scp,1),:) ;
S2 = fftshift(fft(Srx),1);
bitr = cat(1,real(S2),imag(S2));
bitR1 = bitr(1:size(bitr)/2,:)>0;
bitR2 = abs(bitr(1:size(bitr)/2,:))>2;
bitR3 = bitr(size(bitr)/2+1:size(bitr,1),:)>0;
bitR4 = abs(bitr(size(bitr)/2+1:size(bitr,1),:))>2;
bitRx = zeros(N,size(snr,2));
bitRx(1:4:N,:) = bitR1;
bitRx(2:4:N,:) = bitR2;
bitRx(3:4:N,:) = bitR3;
bitRx(4:4:N,:) = bitR4;

% time and frequancy steps
f=((-N/2:1:(N)/2-1)*F).';
fs=((0:N-1)*F).';
t=((0:Ns-1)*T).';
tcp=((0:Ns+Ncp-1)*T).';


% compare rx and tx
cmp = bits==bitRx(1:N,:);
BER = 1-mean(cmp,1);


%plots
figure (1)
 for i=1:size(snr,2)
    clf;
    subplot(2,1,1)
    hold on
    scatter(real(S2(:,i)),imag(S2(:,i)),'filled');
    grid on
    grid minor
    xlabel('I');
    ylabel('Q');
    title('IQ with SNR',snr(i));
    xlim([-4 4])
    ylim([-4 4])
    hold off
    


    pause(0.05)
 end

figure (2)
    semilogy(snr,BER);
    grid on
    grid minor
    xlabel('SNR,db');
    ylabel('BER');




