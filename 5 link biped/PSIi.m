%% Calling this function gives out Psi 
%%=PSIi(x1,x2,alphaa,epsilon)
function p=PSIi(x1,x2,alphaa,epsilon)
x2=epsilon*x2;
phi=x1+(1/(2-alphaa)).*(sign(x2)).*(abs(x2)).^(2-alphaa);
p=-(sign(x2)).*(abs(x2)).^alphaa -(sign(phi)).*(abs(phi)).^(alphaa/(2-alphaa));
p=(1/epsilon^2).*p;
end