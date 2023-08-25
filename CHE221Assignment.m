
% temperature scale %
T = [159.40 153.80 149.40 142.20 133.80 128.30 126.30 122.20 120.30 120.00 119.70 112.70 120.00] ;
% mole fraction of component 1 %
x1i = [0.0870 0.1180 0.1250 0.2190 0.2750 0.4080 0.4800 0.5900 0.6450 0.6510 0.7400 0.8010 0.8840 ];
% mole fraction of compound 2
x2i = 1-x1i;
%vapour mole fraction of component 1
y1exp = [0.5120 0.6210 0.6250 0.7850 0.8070 0.8720 0.8910 0.9160 0.9280 0.9260 0.9460 0.9540 0.9750];
% size of  temperature scale %
[m,n] = size(T);
y1=zeros(m,n);
y2=zeros(m,n);
% plot(x1,y1);
% hold on;
% plot(1-x1,1-y1);
% hold off;
%    data given   %
P = 740.00;
A12 = 123.2046;
A21 = 624.0752;
gamma1inf = 2.12;
gamma2inf = 2.46;
V1L = 106.85;
V2L = 83.14;
R = 8.314;
i=1;
while(i<=n)
    lamma12 = V2L/V1L*(exp(-A12/R*T(1,i)));
    lamma21 = V1L/V2L*(exp(-A21/R*T(1,i)));
    E = -log(x1i(1,i)+(1-x1i(1,i))*lamma12);
    F = -log((1-x1i(1,i))+lamma21*x1i(1,i));
    G = (1-x1i(1,i))*(lamma12/(x1i(1,i)+lamma12*(1-x1i(1,i))) - lamma21/((1-x1i(1,i)) + lamma21*x1i(1,i)));
    H = -x1i(1,i)*(lamma12/(x1i(1,i)+lamma12*(1-x1i(1,i))) - lamma21/((1-x1i(1,i)) + lamma21*x1i(1,i)));
    gamma1 = exp(E+G);
    gamma2 = exp(F+H);

%Antoine constant
    %for toulene
        At = 6.99049;
        Bt = 1346.860;
        Ct = 217.101;
    P1sat = 10^(At-Bt/(T(1,i)+Ct));
    y1(1,i) = x1i(1,i)*gamma1*P1sat/P ;
    
    %for phenol
        Ap = 6.93051;
        Bp = 1382.650;
        Cp = 159.493;
    P2sat = 10^(Ap-Bp/(T(1,i)+Cp));
    y2(1,i) = x2i(1,i)*gamma2*P2sat/P ;

    i=i+1;

end
plot(x1i,T);
hold on;
plot(y1,T);
hold off
