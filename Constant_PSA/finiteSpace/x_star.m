function f = x_star(beta1,y1,beta2)
f = beta1*y1/(1+beta1*y1+beta2*(1-y1));
end