function z = interp2_lim(X,Y,Z,x,y) %implement an input-limited version of interp2
x = min(x,max(X));
x = max(x,min(X));
y = min(y,max(Y));
y = max(y,min(Y));
z = interp2(X,Y,Z,x,y);
end