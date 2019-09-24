%% Refer Asymptotically stable walking for biped robots: analysis via systems with impulse effects
% https://ieeexplore.ieee.org/document/898695
%By Abhijeet Kulkarni amkulk@udel.edu
%%
clear
[m,MH,MT,r,L,g,theta3d,theta1d,alphaa,epsilon,dtheta1_]=simulationparameters();
x0=sigmaa(theta1d,theta3d,dtheta1_);   
x0=impacttransition(x0);
t_span=[0 10];
mechanicalmodeljacoexp;             %calculate expression of impact model
options = odeset('Events',@events);
z1=0;                               %initial foot place
movieVector=[];                     
qall=[];                            %storing q values over simulation
tall=[0];                            %storing t values over simulation
nos=1;
for s=1:4       %number of steps
    [t,x,te,ye,ie]=ode45(@callmechanicalfun,t_span,x0,options);         %call the mechanical function indirectly as ode45 cannot give individual q
    qall=[qall;x];
    tall=[tall;tall(nos)+t];
    p1=zeros(length(t),2);
    p1(:,1)=z1.*ones(length(t),1);
    p2=p1+[(r*(sin(x(:,1))+sin(-x(:,2)))) (r*(cos(x(:,1))-cos(-x(:,2))))];
    p3=p1+[(r*sin(x(:,1))) (r*cos(x(:,1)))];
    p4=p3+[(L*sin(x(:,3))) (L*cos(x(:,3)))];
    z1=p2(length(t),1);
%     x(length(t),:)
    x0=impacttransition(x(length(t),:));
    
    figh=figure(1);
%     length(t)
    for i=1:length(t)
        clf
        hold on;
        plot([p1(i,1) p3(i,1)], [p1(i,2) p3(i,2)] );
        plot([p3(i,1) p2(i,1)], [p3(i,2) p2(i,2)] );
        plot([p3(i,1) p4(i,1)], [p3(i,2) p4(i,2)] );
        xlim([-1 6]);
        ylim([-2 2]);
        grid on;
        title(['t=',num2str(t(i)),'sec']);
        movieVector =[movieVector; getframe(figh)];
    end
    nos=nos+length(t);
end
fig2=figure('Name','thetas');
plot(tall(2:end),qall(:,1),tall(2:end),qall(:,2),tall(2:end),qall(:,3));
legend('theta1','theta2','theta3');
xlabel('time');
ylabel('thetas');
fig3=figure('Name','omegas');
plot(tall(2:end),qall(:,4),tall(2:end),qall(:,5),tall(2:end),qall(:,6));
legend('w1','w2','w3');
xlabel('time');
ylabel('w')
fig4=figure('Name','map');
plot3(qall(:,1),qall(:,4),qall(:,6));
legend;
xlabel('theta1');
ylabel('w1');
zlabel('w3');
fig5=figure('Name','height of swing leg');
plot(tall(2:end),(r.*(cos(qall(:,1))-cos(qall(:,2)))));
xlabel('time');
ylabel('swing foot height');
%% Save Movie
%Create a VideoWriter object and set properties
% myWriter = VideoWriter('moviebip1');            %create an .avi file
myWriter = VideoWriter('moviemp4bip','MPEG-4');   %create an .mp4 file
myWriter.FrameRate = 60;

%Open the VideoWriter object, write the movie, and close the file
open(myWriter);
writeVideo(myWriter, movieVector);
close(myWriter);

function [value,isterminal,direction] = events(t,q)
% Locate the time when height passes through zero in a decreasing direction
% and stop integration.
% 
% value = (1*(cos(q(1))-cos(q(2))));     % detect height = 0
% isterminal = 1;   % stop the integration
% direction = 0;   % negative direction
[m,MH,MT,r,L,g,theta3d,theta1d,alphaa,epsilon,dtheta1_]=simulationparameters();
% if (theta1d-q(1))==0 || q(1)>pi/2 || q(1)<-pi/2
%     value=0;
%     isterminal = 1;
%     direction = 0;
% else
%     value=1;
% end
value=theta1d-q(1);
isterminal=1;
direction=0;
end
