lT=0.625;
lf=0.4;
lt=0.4;
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

fix = open('fixed_point_mod.mat');
ival = fix.x;
x0 = ival(1:10)';
[a,thm,thp] = extractval(ival')
options = odeset('Events',@events,'RelTol',1e-5,'AbsTol',1e-4);
% x0=deg2rad([180 180 0  0 10 0 0 0 0 0]')
t=1;
    z1=0;z2=0;
for step=1:5
     [t,y,te,ye,ie]=ode45(@(t,x)mechfun(t,x,a,thm,thp),[0:0.01:10],x0,options);
      x0=impacttransition(ye);
    % y=x0';
    figh=figure(2);
    size(t)

    for i=1:length(t)
        
        clf
        hold on;
        plot([z1 p_4x(y(i,1),y(i,2),y(i,3),y(i,4),y(i,5),z1,z2)],[z2 p_4y(y(i,1),y(i,2),y(i,3),y(i,4),y(i,5),z1,z2)]);
        plot([p_4x(y(i,1),y(i,2),y(i,3),y(i,4),y(i,5),z1,z2) p_2x(y(i,1),y(i,2),y(i,3),y(i,4),y(i,5),z1,z2)],[p_4y(y(i,1),y(i,2),y(i,3),y(i,4),y(i,5),z1,z2) p_2y(y(i,1),y(i,2),y(i,3),y(i,4),y(i,5),z1,z2)]);
        plot([p_2x(y(i,1),y(i,2),y(i,3),y(i,4),y(i,5),z1,z2) p_1x(y(i,1),y(i,2),y(i,3),y(i,4),y(i,5),z1,z2)],[p_2y(y(i,1),y(i,2),y(i,3),y(i,4),y(i,5),z1,z2) p_1y(y(i,1),y(i,2),y(i,3),y(i,4),y(i,5),z1,z2)]);
        plot([p_2x(y(i,1),y(i,2),y(i,3),y(i,4),y(i,5),z1,z2) p_3x(y(i,1),y(i,2),y(i,3),y(i,4),y(i,5),z1,z2)],[p_2y(y(i,1),y(i,2),y(i,3),y(i,4),y(i,5),z1,z2) p_3y(y(i,1),y(i,2),y(i,3),y(i,4),y(i,5),z1,z2)]);
        plot([p_3x(y(i,1),y(i,2),y(i,3),y(i,4),y(i,5),z1,z2) p_5x(y(i,1),y(i,2),y(i,3),y(i,4),y(i,5),z1,z2)],[p_3y(y(i,1),y(i,2),y(i,3),y(i,4),y(i,5),z1,z2) p_5y(y(i,1),y(i,2),y(i,3),y(i,4),y(i,5),z1,z2)]);
        xlim([-2 2]);
        ylim([-2 2]);
        title(['At time t=',num2str(t(i)),'sec']);
          drawnow;
%         movieVector(j) = getframe(figh);
    end
    z1=p_5x(ye(1),ye(2),ye(3),ye(4),ye(5),z1,z2);
end

% %Create a VideoWriter object and set properties
% myWriter = VideoWriter('moviedoublepend');            %create an .avi file
% % myWriter = VideoWriter('curve','MPEG-4');   %create an .mp4 file
% myWriter.FrameRate = 20;
% 
% %Open the VideoWriter object, write the movie, and close the file
% open(myWriter);
% writeVideo(myWriter, movieVector);
% close(myWriter);


function [value,isterminal,direction] = events(t,q)
% Locate the time when height passes through zero in a decreasing direction
% and stop integration.
lT=0.625;
lf=0.4;
lt=0.4;
value=lt*cos(pi-q(1)-q(2)-q(4))+lf*cos(pi-q(1)-q(2))-lf*cos(q(1)+q(3)-pi)-lt*cos(q(1)+q(3)+q(5)-pi);
isterminal=1;
direction=0;
end

