clear
global g lT lf lt MT Mf Mt IT If It pMT pMf pMt thetap thetam a
g=9.8;
lT=0.625;
lf=0.4;
lt=0.4;
MT = 20;
Mf = 6.8;
Mt = 3.2;
IT = 2.22;
If = 1.08;
It = 0.93;
pMT = 0.2;
pMf = 0.163;
pMt = 0.128;
syms q1 q2 q3 q4 q5 z1 z2 
p_4x(q1,q2,q3,q4,q5,z1,z2) = z1-(lt)*sin(pi-q1-q2-q4);
p_2x(q1,q2,q3,q4,q5,z1,z2)=z1-lt*sin(pi-q1-q2-q4)-(lf)*sin(pi-q1-q2);
p_1x(q1,q2,q3,q4,q5,z1,z2)=z1-lt*sin(pi-q1-q2-q4)-lf*sin(pi-q1-q2)+lT*sin(q1);
p_3x(q1,q2,q3,q4,q5,z1,z2)=z1-lt*sin(pi-q1-q2-q4)-lf*sin(pi-q1-q2)-lf*sin(q1+q3-pi);
p_5x(q1,q2,q3,q4,q5,z1,z2)=z1-lt*sin(pi-q1-q2-q4)-lf*sin(pi-q1-q2)-lf*sin(q1+q3-pi)-lt*sin(q1+q3+q5-pi);
p_4y(q1,q2,q3,q4,q5,z1,z2) = z2+(lt)*cos(pi-q1-q2-q4);
p_2y(q1,q2,q3,q4,q5,z1,z2)=z2+lt*cos(pi-q1-q2-q4)+(lf)*cos(pi-q1-q2);
p_1y(q1,q2,q3,q4,q5,z1,z2)=z2+lt*cos(pi-q1-q2-q4)+lf*cos(pi-q1-q2)+lT*cos(q1);
p_3y(q1,q2,q3,q4,q5,z1,z2)=z2+lt*cos(pi-q1-q2-q4)+lf*cos(pi-q1-q2)-lf*cos(q1+q3-pi);
p_5y(q1,q2,q3,q4,q5,z1,z2)=z2+lt*cos(pi-q1-q2-q4)+lf*cos(pi-q1-q2)-lf*cos(q1+q3-pi)-lt*cos(q1+q3+q5-pi);

% fix = open('new_value_mod.mat');
fix = open('good_fixed_point.mat');
% fix = open('fixed_point_mod.mat');
ival = fix.x;
x0 = ival(1:10)';
[a,thetap,thetam] = extractval(ival');
options = odeset('Events',@eventss,'RelTol',1e-10,'AbsTol',1e-10);
% x0=deg2rad([180 180 0  0 10 0 0 0 0 0]')
t=0;
te=0;
k=0;
z1=-1.5;z2=0;
yout=[];
tout=[];
figh=figure('Name','Animation');
for step=1:10
     [t,y,te,ye,ie]=ode45(@(t,x)mechfun(t,x),[te:0.01:10],x0,options);
      x0=impacttransition(ye);
    % y=x0';
    step
    tout=[tout;t];
    yout=[yout;y];
if false
    for i=1:length(t)
        k=k+1;
        clf
        hold on;
        plot([z1 p_4x(y(i,1),y(i,2),y(i,3),y(i,4),y(i,5),z1,z2)],[z2 p_4y(y(i,1),y(i,2),y(i,3),y(i,4),y(i,5),z1,z2)]);
        plot([p_4x(y(i,1),y(i,2),y(i,3),y(i,4),y(i,5),z1,z2) p_2x(y(i,1),y(i,2),y(i,3),y(i,4),y(i,5),z1,z2)],[p_4y(y(i,1),y(i,2),y(i,3),y(i,4),y(i,5),z1,z2) p_2y(y(i,1),y(i,2),y(i,3),y(i,4),y(i,5),z1,z2)]);
        plot([p_2x(y(i,1),y(i,2),y(i,3),y(i,4),y(i,5),z1,z2) p_1x(y(i,1),y(i,2),y(i,3),y(i,4),y(i,5),z1,z2)],[p_2y(y(i,1),y(i,2),y(i,3),y(i,4),y(i,5),z1,z2) p_1y(y(i,1),y(i,2),y(i,3),y(i,4),y(i,5),z1,z2)]);
        plot([p_2x(y(i,1),y(i,2),y(i,3),y(i,4),y(i,5),z1,z2) p_3x(y(i,1),y(i,2),y(i,3),y(i,4),y(i,5),z1,z2)],[p_2y(y(i,1),y(i,2),y(i,3),y(i,4),y(i,5),z1,z2) p_3y(y(i,1),y(i,2),y(i,3),y(i,4),y(i,5),z1,z2)]);
        plot([p_3x(y(i,1),y(i,2),y(i,3),y(i,4),y(i,5),z1,z2) p_5x(y(i,1),y(i,2),y(i,3),y(i,4),y(i,5),z1,z2)],[p_3y(y(i,1),y(i,2),y(i,3),y(i,4),y(i,5),z1,z2) p_5y(y(i,1),y(i,2),y(i,3),y(i,4),y(i,5),z1,z2)]);
        xlim([-2 2]);
        ylim([-0.2 2]);
        pbaspect([1 1 1])
        axis on;
        grid on;
        legend('StanceLeg','StanceFemur','Torso','SwingFemur','SwingLeg')
        title(['At time t=',num2str(t(i)),'sec']);
        drawnow;
          
          
        movieVector(k) = getframe(figh);
    end
end
    z1=p_5x(ye(1),ye(2),ye(3),ye(4),ye(5),z1,z2);
end
% yo
x(1:10) = x0;
x(11:26) = ival(11:26);
save('new_fixed_point_mod.mat','x');

% for i=1:length(tout)
% [Q,u(i)]=mechfun(tout(i),yout(1,:)',a,thm,thp);
% fcon=contact_force(yout,u);
% end
% plot(tout,fcon(:,2));
th = yout(:,1) + yout(:,2) + 0.5*yout(:,4);
thd = yout(:,6) + yout(:,7) + 0.5*yout(:,9);
limitcycle=figure('Name','Limit cycle');
plot(th,thd)
saveas(limitcycle,'limitcycle.jpg');
% [t,y,te,ye,ie]=ode45(@(t,x)mechfun(t,x,a,thm,thp),[0 10],ival(1:10),options);
t=tout;
y=yout;
for j = 1:length(t)
            [Q,u(j,:)] = mechfun(t,y(j,:)',a,thm,thp);
            Fcontact(j,:) = contact_force(y(j,:)',u(j,:)');
             ftfn(j)=Fcontact(j,1)/Fcontact(j,2);
end
FtFn=figure('Name','Ft/Fn');
hold on;
plot(t,Fcontact(:,1));
plot(t,Fcontact(:,2));
legend('1','2')
mu=figure('Name','Ft/Fn')
plot(t,Fcontact(:,1)./Fcontact(:,2))
hold off;
saveas(mu,'mu.jpg');
saveas(FtFn,'FtFn.jpg');
intorque=figure('Name','Input torques');
hold on;
plot(t,u(:,1))
plot(t,u(:,2))
plot(t,u(:,3))
plot(t,u(:,4))
legend('u1','u2','u3','u4')
hold off;
saveas(intorque,'inputtorque.jpg');
if false
%Create a VideoWriter object and set properties
myWriter = VideoWriter('moviedoublepend');            %create an .avi file
% myWriter = VideoWriter('curve','MPEG-4');   %create an .mp4 file
myWriter.FrameRate = 20;

%Open the VideoWriter object, write the movie, and close the file
open(myWriter);
writeVideo(myWriter, movieVector);
close(myWriter);
end


% function [value,isterminal,direction] = events(t,q)
% % Locate the time when height passes through zero in a decreasing direction
% % and stop integration.
% lT=0.625;
% lf=0.4;
% lt=0.4;
% value=lt*cos(pi-q(1)-q(2)-q(4))+lf*cos(pi-q(1)-q(2))-lf*cos(q(1)+q(3)-pi)-lt*cos(q(1)+q(3)+q(5)-pi);
% isterminal=1;
% direction=0;
% end

