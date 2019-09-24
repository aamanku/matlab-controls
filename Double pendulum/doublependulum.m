clear
clc
close all
% equations from https://www.astro.umd.edu/~adhabal/V1/Reports/Order_and_Chaos.pdf
% refer for animation - http://faculty.washington.edu/lum/EducationalVideoFiles/Matlab05/AnimationMatlab.m
% also refer for decoupling https://diego.assencio.com/?index=1500c66ae7ab27bb0106467c68feebc6
%% Initial
l1=1;
l2=1;
x0=[2;0;1;0];
[t,x]=ode45(@fun,[0 10],x0);


%% Display part   
hold on;
path=figure(1);
plot(l1.*sin(x(:,1)),-l1.*cos(x(:,1)));
plot(l1.*sin(x(:,1))+l2.*sin(x(:,3)),-l1.*cos(x(:,1))-l2.*cos(x(:,3)));
saveas(path,'pathdoublependulum.jpg');

figh=figure(2);

for i=1:length(t)
% plot(cos(x(i,1)),sin(x(i,1)))
    clf
    hold on; 
    plot([0 l1*sin(x(i,1))],[0 -l1*cos(x(i,1))]);
    plot([l1*sin(x(i,1)) l1*sin(x(i,1))+l2*sin(x(i,3))],[-cos(x(i,1)) -cos(x(i,1))-l2*cos(x(i,3))]);
    plot(l1*sin(x(i,1)),-l1*cos(x(i,1)),'o');
    plot(l1*sin(x(i,1))+l2*sin(x(i,3)),-l1.*cos(x(i,1))-l2*cos(x(i,3)),'o');
    grid on;
    xlabel('x');
    ylabel('y');
    xlim([-2 2]);
    ylim([-2 2]);
    title(['At time t=',num2str(t(i)),'sec']);
%     drawnow;
    movieVector(i) = getframe(figh);
end

%% Save Movie
%Create a VideoWriter object and set properties
myWriter = VideoWriter('moviedoublepend');            %create an .avi file
% myWriter = VideoWriter('curve','MPEG-4');   %create an .mp4 file
myWriter.FrameRate = 20;

%Open the VideoWriter object, write the movie, and close the file
open(myWriter);
writeVideo(myWriter, movieVector);
close(myWriter);
 
%% Functions
function f=fun(t,x)
c1=2.16;
c2=0.27;
c3=1;
c4=172.5;
c5=32.6;
f=[x(2);(2*c2*c4*sin(x(1)) + (c3^2)*(x(2)^2)*sin(x(1)-x(3))*cos(x(1)-x(3))+2*c2*c3*(x(4)^2)*sin(x(1)-x(3)) - c3*c5*sin(x(3))*cos(x(1)-x(3))) / ((c3^2)*cos(x(1)-x(3))^2-4*c1*c2);x(4);(2*c1*c5*sin(x(3)) - (c3^2)*(x(4)^2)*sin(x(1)-x(3))*cos(x(1)-x(3))-2*c1*c3*(x(2)^2)*sin(x(1)-x(3)) - c3*c4*sin(x(1))*cos(x(1)-x(3))) / ((c3^2)*cos(x(1)-x(3))^2-4*c1*c2)];

% p1dot=(2*c2*c4*sin(x(1)) + (c3^3)*(x(2)^2)*sin(x(1)-x(3))*cos(x(1)-x(3))
%  2*c2*c3*(x(4)^2)*sin(x(1)-x(3)) - c3*c5*sin(x(3))*cos(x(1)-x(3))) / ((c3^2)*cos(x(1)-x(3))^2-4*c1*c2)
% 
% p2dot=(2*c2*c5*sin(x(3)) + (c3^3)*(x(4)^2)*sin(x(1)-x(3))*cos(x(1)-x(3))
%  2*c1*c3*(x(2)^2)*sin(x(1)-x(3)) - c3*c4*sin(x(1))*cos(x(1)-x(3))) / ((c3^2)*cos(x(1)-x(3))^2-4*c1*c2)
end