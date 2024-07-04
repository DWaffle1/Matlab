N=4096;
F=1000;
T=1/F;


% creating tx signal
bits = rand(N,1)>0.5;
Ibits = 2*bits(1:2:N)-1;
Qbits = 2*bits(2:2:N)-1;
S = (Ibits+1i*Qbits)*1/sqrt(2);
Ns=size(S,1);
Ncp=round(Ns/10);
Stx = ifft(fftshift(S));
Scp = Stx(Ns-Ncp+1:Ns);
NewS = cat(1,Scp,Stx);


% noise
snr = [41:-0.5:-10];
N0 = 10.^(-snr./10);
noise = sqrt((N0)/2).*randn(Ns+Ncp,1)+1i*sqrt((N0)/2).*randn(Ns+Ncp,1);

NewS = NewS+noise;

% processing rx signal
Srx = NewS(size(Scp,1)+1:size(Stx,1)+size(Scp,1),:) ;
S2 = fftshift(fft(Srx),1);
b=zeros(N,size(snr,2));
b(1:2:N,:)=real(S2);
b(2:2:N,:)=imag(S2);
bitsRx=b>0;


% time and frequancy steps
f=((-N/2:1:(N)/2-1)*F).';
fs=((0:N-1)*F).';
t=((0:Ns-1)*T).';
tcp=((0:Ns+Ncp-1)*T).';

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
