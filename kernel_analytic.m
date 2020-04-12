function [F]=kernel_analytic(sst,XY)
% Analytic kernel from Del Pezzo et al. 2016
deltaxy=0.2;
xx=XY(:,1);
yy=XY(:,2);
D1=sqrt((sst(4)-sst(1))^2+(sst(5)-sst(2))^2);
F1=1/2/pi/deltaxy^2/D1^2;
F2=1/deltaxy^2/D1^2;
F3=F1*(0.5*exp(-abs((xx-(sst(1)+sst(4))/2).^2*F2/2+...
    (yy-(sst(2)+sst(5))/2).^2*F2/0.5))+...
    exp(-abs((xx-sst(1)).^2*F2/2+...
    (yy-sst(2)).^2*F2/2))+...
    exp(-abs((xx-sst(4)).^2*F2/2+...
    (yy-sst(5)).^2*F2/2)));
if find(F3)>0
    F=F3/sum(F3);
else
    F=F3;
end
no=F<0.0001;
F(no)=0;
