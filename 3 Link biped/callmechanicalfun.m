%% This function calls the mechanicalfun file with proper inputs
function y=callmechanicalfun(~,q)
y=mechanicalfun(q(1),q(2),q(3),q(4),q(5),q(6));
end