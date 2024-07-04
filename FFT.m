F=1000;
T=1/F;
N=100;
t=([0:1:N-1].')*T;
f=([-N/2:1:(N-1)/2].')*F/N;

S1 =[ zeros(floor(N/3),1 ); ones(N-2*floor(N/3),1); zeros(floor(N/3),1 )];
FT1=fft(S1);

S2 =[ zeros(floor(4*N/9),1 ); ones(N-2*floor(4*N/9),1); zeros(floor(4*N/9),1 )];
FT2=fft(S2);

S3 =[ zeros(floor(N/9),1 ); ones(N-2*floor(N/9),1); zeros(floor(N/9),1 )];
FT3=fft(S3);

figure
    subplot(3,2,1)
    hold on
    stem(t,S1);
    grid on
    grid minor
    xlabel('Время,с');
    ylabel('Сигнал');
    title('Сигнал от времени 1(средний)');
    hold off

    subplot(3,2,2)
    hold on
    stem(f,fftshift(abs(FT1))/N)
    grid on
    grid minor
    xlabel('Частота, Гц');
    ylabel('Сигнал');
    title('Спектр сигнала 1(средний)');
    hold off

    subplot(3,2,3)
    hold on
    stem(t,S2,'red');
    grid on
    grid minor
    xlabel('Время,с');
    ylabel('Сигнал');
    title('Сигнал от времени 2(узкий)');
    hold off

    subplot(3,2,4)
    hold on
    stem(f,fftshift(abs(FT2))/N,'red')
    grid on
    grid minor
    xlabel('Частота, Гц');
    ylabel('Сигнал');
    title('Спектр сигнала 2(узкий)');
    hold off

    subplot(3,2,5)
    hold on
    stem(t,S3,'green');
    grid on
    grid minor
    xlabel('Время,с');
    ylabel('Сигнал');
    title('Сигнал от времени 3(широкий)');
    hold off

    subplot(3,2,6)
    hold on
    stem(f,fftshift(abs(FT3))/N,'green')
    grid on
    grid minor
    xlabel('Частота, Гц');
    ylabel('Сигнал');
    title('Спектр сигнала 3(широкий)');
    hold off
