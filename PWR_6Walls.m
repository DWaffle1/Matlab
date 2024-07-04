%%%%%%%%%%%%%% set initial conditions

C_Rsrc_S = rand(3,1);
C_Rdst_S = rand(3,1);
g=100000;
time = 1:g;
Coord_Rsrc =[1;1;1];                %[ (C_Rsrc_S(1)*5)-2.5 ; (C_Rsrc_S(2)*6)-3 ; C_Rsrc_S(3)*3];
Coord_Rdst =Coord_Rsrc+[0.1;0;0]*time;
Coord_Rdst=reshape(Coord_Rdst,[3,1,g]);


% [ (C_Rdst_S(1)*5)-2.5 ; (C_Rdst_S(2)*6)-3 ; C_Rdst_S(3)*3]
    % mask=[-2.5;-3;0]<Coord_Rdst & Coord_Rdst<[2.5;3;3];
    % Coord_Rdst=Coord_Rdst.*mask;

%%%%%%%%%%%%%%

%%%%%%%%%%%%%% create a structure for describing walls
stenka=struct("name","floor","t1",[1;0;0],"t2",[0;1;0],"PlaneDot",[-2.5 ; -3 ; 0],"Size1",5,"Size2",6,"Norm",[0;0;1]);
stenka(end+1).name = "wallYD";
stenka(end).t1 = [1;0;0];
stenka(end).t2 = [0;0;1];
stenka(end).PlaneDot = [-2.5 ; -3 ; 0];
stenka(end).Size1 = 5;
stenka(end).Size2 = 3;
stenka(end).Norm = [0;1;0];

stenka(end+1).name = "wallXL";
stenka(end).t1 = [0;0;1];
stenka(end).t2 = [0;1;0];
stenka(end).PlaneDot = [-2.5 ; -3 ; 0];
stenka(end).Size1 = 3;
stenka(end).Size2 = 6;
stenka(end).Norm = [1;0;0];

stenka(end+1).name = "roof";
stenka(end).t1 = [1;0;0];
stenka(end).t2 = [0;1;0];
stenka(end).PlaneDot = [-2.5 ; -3 ; 3];
stenka(end).Size1 = 5;
stenka(end).Size2 = 6;
stenka(end).Norm = [0;0;1];

stenka(end+1).name = "wallYU";
stenka(end).t1 = [1;0;0];
stenka(end).t2 = [0;0;1];
stenka(end).PlaneDot = [-2.5 ; 3 ; 0];
stenka(end).Size1 = 5;
stenka(end).Size2 = 3;
stenka(end).Norm = [0;1;0];

stenka(end+1).name = "wallXR";
stenka(end).t1 = [0;0;1];
stenka(end).t2 = [0;1;0];
stenka(end).PlaneDot = [2.5 ; -3 ; 0];
stenka(end).Size1 = 3;
stenka(end).Size2 = 6;
stenka(end).Norm = [1;0;0];

t1=cat(2,stenka(:).t1);
t2=cat(2,stenka(:).t2);
CP=cat(2,stenka(:).PlaneDot);
size1=cat(1,stenka(:).Size1);
size2=cat(1,stenka(:).Size2);
n=cat(2,stenka(:).Norm);
%%%%%%%%%%%%%%

%%%%%%%%%%%%%% initialization of functions

Coord_Rint = intersection(Coord_Rdst,Coord_Rsrc,CP,n);
D=plosk(6 ,CP,size1,size2,t1,t2);
Coord_Ridst6=ri(Coord_Rdst,CP,n); % imaginary Rdst for 6 walls
Coord_Ridst=ri(Coord_Rdst,CP(:,1),n(:,1));  % imaginary Rdst for 1 wall



%%%%%%%%%%%%%% Calculating power

Rsid=(sqrt(sum( (Coord_Ridst - Coord_Rsrc) .^2,1)));    % distance between Src and imaginary Dst
Rsid6=(sqrt(sum( (Coord_Ridst6 - Coord_Rsrc) .^2,1))); % distance in case of 6 walls
Rsd=sqrt(sum((Coord_Rdst-Coord_Rsrc).^2,1));    % distance between Src and real Dst

SD=exp(1i*2*pi.*Rsd/0.01)./Rsd;      
SID=exp((1i*2*pi.*Rsid/0.01)+1i*pi )./Rsid;
SID6=exp((1i*2*pi.*Rsid6/0.01)+1i*pi )./Rsid6;

%Calculating power for 1 wall
res1=cat(2,SID,SD);
res1=sum(res1,2);
res1=abs(res1).^2;
res1=reshape(res1,[g,1]);

%Calculating power for 6 walls

res6=cat(2,SID6,SD);
res6=sum(res6,2);
res6=abs(res6).^2;
res6=reshape(res6,[g,1]);

%Calculating power without walls
resSD=abs(SD).^2;
resSD=reshape(resSD,[g,1]);

%%%%%%%%%%%%%% drawing plot
% figure
%  for j=1:g 
%     clf;
%     hold on
%     view(55,16)
%     lineSD = [Coord_Rsrc Coord_Rdst(:,j)];
%     plot3(lineSD(1,:) , lineSD(2,:) , lineSD(3,:) );
%     plot3(Coord_Rsrc(1),Coord_Rsrc(2),Coord_Rsrc(3),'ro');
%     plot3(Coord_Rdst(1,j),Coord_Rdst(2,j),Coord_Rdst(3,j),'bo');
%     for i=1:1:6
%         fill3(D(1,:,i),D(2,:,i),D(3,:,i), [1,1,1],'FaceAlpha',0.1)
%         lineSI = [Coord_Rsrc Coord_Rint(:,i,j)];
%         lineID = [Coord_Rint(:,i,j) Coord_Rdst(:,j)];
%         plot3(lineSI(1,:) , lineSI(2,:) , lineSI(3,:) );
%         plot3(lineID(1,:) , lineID(2,:) , lineID(3,:) );
%         plot3(Coord_Rint(1,i,j),Coord_Rint(2,i,j),Coord_Rint(3,i,j),'go');
%     end
%     pause(0.0001);
%     axis equal    
%     hold off 
%  end


%%%%%%%%%%%%%% drawing plot for power

figure
 loglog(time,res1)
 hold on
 loglog(time,resSD)

figure
 loglog(time,res6)
 hold on
 loglog(time,resSD)

%%%%%%%%%%%%%%

function D = plosk(count ,CP,size1,size2,t1,t2) % counting coordinates of wall's corners  
    for i=1:count
        D(:,1,i) = CP(:,i);
        D(:,2,i) = CP(:,i) + (size1(i) .* t1(:,i));
        D(:,3,i) = CP(:,i) + (size1(i) .* t1(:,i)) + (size2(i) .* t2(:,i));
        D(:,4,i) = CP(:,i) + (size2(i) .* t2(:,i)); 
    end
end


function Coord_Rint = intersection(Coord_Rdst,Coord_Rsrc,CP,n) % finding point of intersection
    Coord_RiDst = Coord_Rdst - 2*sum((Coord_Rdst - CP) .* n, 1).*n;
    V = (Coord_RiDst - Coord_Rsrc) ./ (sqrt(sum( (Coord_RiDst - Coord_Rsrc) .^2,1)));
    T = - (sum(  (Coord_Rsrc - CP).*n , 1)) ./ (sum ( n .* V , 1));
    Coord_Rint = (V .* T) + Coord_Rsrc;
end


function Coord_RiDst = ri(Coord_Rdst,CP,n) % finding point of Ridst
    Coord_RiDst = Coord_Rdst - 2*sum((Coord_Rdst - CP) .* n, 1).*n;
end



