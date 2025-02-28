function y = interp1_lim(X,Y,x) %implement an input-limited version of interp1
x = min(x,max(X));
x = max(x,min(X));
y = interp1(X,Y,x);
end