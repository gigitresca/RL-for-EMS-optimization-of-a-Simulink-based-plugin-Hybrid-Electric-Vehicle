function z = interp3_lim(X,Y,W,Z,x,y,w) %implement an input-limited version of interp2
x = min(x,max(X));
x = max(x,min(X));
y = min(y,max(Y));
y = max(y,min(Y));
w = min(w,max(W));
w = max(w,min(W));
z = interp3(X,Y,W,Z,x,y,w);
end