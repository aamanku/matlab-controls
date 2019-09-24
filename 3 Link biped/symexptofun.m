syms x y z t
r = (x + y/2 + z/3)*exp(-t);
matlabFunction(r,'File','myfile','Vars',{t,[x y z]});