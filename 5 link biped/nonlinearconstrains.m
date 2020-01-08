function [C,Ceq]=nonlinearconstrains(qpoint)
qpoint=qpoint(:);
%qpoint-till 10 states \11 to 26 coeficients of bezeir
MAX_TORQUE=170;
MAX_FRICTIONCOEFFICIENT=1;
MAX_AVERAGESPEED=10;
MIN_AVERAGESPEED=0.7;
MIN_STEPLENGTH=0.2;

[a,thm,thp] = extractval(qpoint);
qpoint=qpoint(:);
options = odeset('Events',@eventss,'RelTol',1e-5,'AbsTol',1e-4);
[t,y,te,ye,ie]=ode45(@(t,x)mechfun(t,x,a,thm,thp),[0 10],qpoint(1:10),options);
if isempty(ye)
    emp=1
    C=999999*ones(1,13).*(1+rand);
    Ceq=999999*(1+rand)*ones(1,10);
%     C=Inf(1,12);
%     Ceq=Inf;
else
    for j = 1:length(t)
            [Q,u(j,:)] = mechfun(t,y(j,:)',a,thm,thp);
            Fcontact(j,:) = conforce(y(j,:)',u(j,:)');
    end
    q4 = y(:,4);
    q5 = y(:,5);
    qp = impacttransition(ye);
    friction_cone = abs(Fcontact(:,1)./Fcontact(:,2));
    step_length = steplength(ye);
    average_speed = step_length/te;
    C=[(abs(qpoint(1))-0.3*pi/2);(qpoint(2)-pi);(qpoint(3)-1.5*pi);(qpoint(4)-pi/2);(qpoint(5)-pi/2)];
    C = [C;max(max(abs(u)))-MAX_TORQUE; 
        -min(q4); 
        -min(q5); 
        max(friction_cone)-MAX_FRICTIONCOEFFICIENT ;
        -average_speed+MIN_AVERAGESPEED ; 
        average_speed-MAX_AVERAGESPEED ; 
        -step_length+MIN_STEPLENGTH;
        thp-thm]';  
    Ceq = (qp - qpoint(1:10))';
end
end
