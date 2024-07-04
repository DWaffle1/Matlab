Coord_Rsrc=[1;1;1];
g=500000;
time=1:1:g;
Coord_Rdst=Coord_Rsrc+[0.1;0;0]*time;

CP=[-2.5 ; -3 ; 0];
n=[0;0;1];


Coord_RiDst = Coord_Rdst - 2*sum((Coord_Rdst - CP) .* n, 1).*n;
V = (Coord_RiDst - Coord_Rsrc) ./ (sqrt(sum( (Coord_RiDst - Coord_Rsrc) .^2,1)));
T = - (sum(  (Coord_Rsrc - CP).*n , 1)) ./ (sum ( n .* V , 1));
Coord_Rint = (V .* T) + Coord_Rsrc;

Rsi=sqrt(sum((Coord_Rint-Coord_Rsrc).^2,1));
Rid=sqrt(sum((Coord_RiDst-Coord_Rint).^2,1));
R=Rsi+Rid;
Rsd=sqrt(sum((Coord_Rdst-Coord_Rsrc).^2,1));
k=(2*pi)/0.01;
SID =exp((1i*k*R)+1i*pi )./(R);
%ID =exp( (1i*k*Rid)+pi )./(Rid);
SD =exp(1i*k*Rsd)./(Rsd);

res=SID+SD;
res=sum(res,1);
res=abs(res).^2;
res=reshape(res,[1,g]);
res2=abs(SD).^2;


figure
    loglog(time,res)
    hold on
    loglog(time,res2)    
