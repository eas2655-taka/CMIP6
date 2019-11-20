function r = nancorr(X,Y);

mask = X+Y;
X(isnan(mask))=[];
Y(isnan(mask))=[];

C=cov(X,Y);
r=C(1,2)/sqrt(C(1,1)*C(2,2));

return;
