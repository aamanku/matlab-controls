function cost=fitnessfunction(q)
q=q(:);
[a,thm,thp] = extractval(q);
options = odeset('Events',@eventss,'RelTol',1e-5,'AbsTol',1e-4);
[t,y,te,ye,ie]=ode45(@(t,x)mechfun(t,x,a,thm,thp),[0 10],q(1:10),options);
if isempty(ye)
    cost=999999*(1+rand);
    efc=1
%     cost=Inf;
else
    for j = 1:length(t)
            [Q,u(j,:)] = mechfun(t,y(j,:)',a,thm,thp);
         
    end
    step_length = steplength(ye);
    cost=abs(abs(norm(u)/step_length));
    
end
end
