N=4096;
Ncp=round(N/10);
F=1000;
T=1/F;


% time and frequancy steps
f=((-N/2:1:(N)/2-1)*F).';
fs=((0:N-1)*F).';
t=((0:N-1)*T).';
tcp=((0:N+Ncp-1)*T).';


% creating tx signal
bits = rand(N,1)>0.5;
S = 2*bits-1;
Stx = ifft(fftshift(S))*sqrt(N);
Scp = Stx(N-Ncp+1:N);
NewS = cat(1,Scp,Stx);


% noise
snr = [30:-0.5:-10];
N0 = 10.^(-snr./10);
noise = sqrt((N0)/2).*randn(N+Ncp,1)+1i*sqrt((N0)/2).*randn(N+Ncp,1);

NewS = NewS+noise;

% processing rx signal
Srx = NewS(size(Scp,1)+1:size(Stx,1)+size(Scp,1),:) ;
S2 = fftshift(fft(Srx)/sqrt(N),1);
bitsRx = real(S2)>0;

% compare rx and tx
cmp = bits==bitsRx(1:N,:);
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
    xlim([-2 2])
    ylim([-2 2])
    hold off
    

    subplot(2,1,2)
    hold on
    semilogy(f,cmp(1:N,i),'o');
    grid on
    grid minor
    xlabel('f,Гц');
    ylabel('Bits');
    title('Comparison');
    hold off
    pause(0.05)
 end

 figure (2)
    semilogy(snr,BER);
    grid on
    grid minor
    xlabel('SNR,db');
    ylabel('BER');


