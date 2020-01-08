warning('off','all');
x = 1+rand(26,1)
LB=x - 4;
UB=x + 3;
options = optimset('Display','iter-detailed','TolFun',1e-20,'TolX',1e-20,'MaxFunEvals',1e7,'UseParallel',true,'PlotFcn',@optimplotfval);
x = fmincon('fitnessfunction',x,[],[],[],[],LB',UB','nonlinearconstrains',options)
save('output.mat','x');