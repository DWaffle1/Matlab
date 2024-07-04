%%%%%%%%%%%%%% set initial conditions

C_Rsrc_S = rand(3,1);
C_Rdst_S = rand(3,1);

Coord_Rsrc = [ (C_Rsrc_S(1)*5)-2.5 ; (C_Rsrc_S(2)*6)-3 ; C_Rsrc_S(3)*3];
Coord_Rdst = [ (C_Rdst_S(1)*5)-2.5 ; (C_Rdst_S(2)*6)-3 ; C_Rdst_S(3)*3];




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

%%%%%%%%%%%%%%

%%%%%%%%%%%%%% drawing plot
hold on
lineSD = [Coord_Rsrc Coord_Rdst];
plot3(lineSD(1,:) , lineSD(2,:) , lineSD(3,:) );
plot3(Coord_Rsrc(1),Coord_Rsrc(2),Coord_Rsrc(3),'ro');
plot3(Coord_Rdst(1),Coord_Rdst(2),Coord_Rdst(3),'bo');
for i=1:1:6
   fill3(D(1,:,i),D(2,:,i),D(3,:,i), [1,1,1],'FaceAlpha',0.1)
   lineSI = [Coord_Rsrc Coord_Rint(:,i)];
   lineID = [Coord_Rint(:,i) Coord_Rdst];
   plot3(lineSI(1,:) , lineSI(2,:) , lineSI(3,:) );
   plot3(lineID(1,:) , lineID(2,:) , lineID(3,:) );
   plot3(Coord_Rint(1,i),Coord_Rint(2,i),Coord_Rint(3,i),'go');
end
axis equal
view(3)
hold off

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
    V = (Coord_RiDst - Coord_Rsrc) ./ (sqrt(sum( (Coord_RiDst - Coord_Rsrc) .^2,1)))
    T = - (sum(  (Coord_Rsrc - CP).*n , 1)) ./ (sum ( n .* V , 1))
    Coord_Rint = (V .* T) + Coord_Rsrc;
end