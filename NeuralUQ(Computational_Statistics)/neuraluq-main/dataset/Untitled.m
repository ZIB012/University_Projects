P = polyfit(X,Y,20);
plot(xx, polyval(P,xx))
hold on
plot(X,Y,'ro')