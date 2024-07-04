N=300;
St=[zeros(N/3,1);ones(N/3,1);zeros(N/3,1)];
n=repmat([0:1:N-1],N,1).* repmat([0:1:N-1].',1,N);

W=( 1/(sqrt(N)) )*exp( (-i*2*pi*n)/N );
Sf=W*St;

t=[-N/2:1:(N/2)-1].';
f=[-N/2:1:(N/2)-1].';

figure
    subplot(2,1,1)
    hold on
    plot(t,St,'b.','MarkerSize',4)
    plot(t,St,'MarkerSize',4)
    grid on
    grid minor
    xlabel('Время');
    ylabel('Сигнал');
    title('Сигнал от времени');
    hold off

    subplot(2,1,2)
    hold on
    plot(f,Sf,'MarkerSize',4)
    grid on
    grid minor
    xlabel('Частота');
    ylabel('Сигнал');
    title('Спектр сигнала');
    hold off


