a=0;b=1;
%dividing this interval in N equal parts
N=16; rangeN=1;
N1=32;rangeN1=2;
N2=64;rangeN2=4;
N3=128;rangeN3=8;
N4=256;rangeN4=16;
N5=512;rangeN5=32;
%h be the steplength
h=(b-a)/N;
h1=(b-a)/N1;
h2=(b-a)/N2;
h3=(b-a)/N3;
h4=(b-a)/N4;
h5=(b-a)/N5;

x=linspace(a,b,N+1);
x1=linspace(a,b,N1+1);
x2=linspace(a,b,N2+1);
x3=linspace(a,b,N3+1);
x4=linspace(a,b,N4+1);
x5=linspace(a,b,N5+1);

%taking the transpose of x
x=x';
x1=x1';
x2=x2';
x3=x3';
x4=x4';
x5=x5';
%r.h.s of the given differential equation
f=sin((x.^2)+1);
f=sin((x1.^2)+1);
f=sin((x2.^2)+1);
f=sin((x3.^2)+1);
f=sin((x4.^2)+1);
f=sin((x5.^2)+1);

%excluding the two boundary condition
X=x(2:N);
X1=x1(2:N1);
X2=x2(2:N2);
X3=x3(2:N3);
X4=x4(2:N4);
X5=x5(2:N5);
%same for the r.h.s
F=sin((X.^2)+1);
F1=sin((X1.^2)+1);
F2=sin((X2.^2)+1);
F3=sin((X3.^2)+1);
F4=sin((X4.^2)+1);
F5=sin((X5.^2)+1);

AX=2.^X;
AX1=2.^X1;
AX2=2.^X2;
AX3=2.^X3;
AX4=2.^X4;
AX5=2.^X5;
%produce the sparse matrix A
n=N-1;e=ones(n,1);
n1=N1-1;e1=ones(n1,1);
n2=N2-1;e2=ones(n2,1);
n3=N3-1;e3=ones(n3,1);
n4=N4-1;e4=ones(n4,1);
n5=N5-1;e5=ones(n5,1);


A=((-1)/(h*h))*spdiags([e -(2+h*h*AX).*e e],-1:1,n,n);
A1=((-1)/(h1*h1))*spdiags([e1 -(2+h1*h1*AX1).*e1 e1],-1:1,n1,n1);
A2=((-1)/(h2*h2))*spdiags([e2 -(2+h2*h2*AX2).*e2 e2],-1:1,n2,n2);
A3=((-1)/(h3*h3))*spdiags([e3 -(2+h3*h3*AX3).*e3 e3],-1:1,n3,n3);
A4=((-1)/(h4*h4))*spdiags([e4 -(2+h4*h4*AX4).*e4 e4],-1:1,n4,n4);
A5=((-1)/(h5*h5))*spdiags([e5 -(2+h5*h5*AX5).*e5 e5],-1:1,n5,n5);


%full(A)
%full(A1)
%%full(A2)
%full(A3)
%full(A4)
%full(A5)

%solving the equation AU=F
U=A\F;
U1=A1\F1;
U2=A2\F2;
U3=A3\F3;
U4=A4\F4;
U5=A5\F5;

format long
%plotting in the diagram
plot(X,U(rangeN:rangeN:end),'o')
hold on
plot(X1(rangeN1:rangeN1:end),U1(rangeN1:rangeN1:end),'*')
plot(X2(rangeN2:rangeN2:end),U2(rangeN2:rangeN2:end),'d')
plot(X3(rangeN3:rangeN3:end),U3(rangeN3:rangeN3:end),'s')
plot(X4(rangeN4:rangeN4:end),U4(rangeN4:rangeN4:end),'^')
plot(X5(rangeN5:rangeN5:end),U5(rangeN5:rangeN5:end),'p')

N, max(U(rangeN:rangeN:end))
N1, max(U1(rangeN1:rangeN1:end)-U(rangeN:rangeN:end))
N2, max(U2(rangeN2:rangeN2:end)-U1(rangeN1:rangeN1:end))
N3, max(U3(rangeN3:rangeN3:end)-U2(rangeN2:rangeN2:end))
N4, max(U4(rangeN4:rangeN4:end)-U3(rangeN3:rangeN3:end))
N5, max(U5(rangeN5:rangeN5:end)-U4(rangeN4:rangeN4:end))

legend('N=16','N=32','N=64','N=128','N=256','N=512')

%% Answers
% N   MaxError
% 16  0.0
% 32  1.789504878942083e-05 
% 64  4.467460969811987e-06 
% 128 1.116472152146164e-06 
% 256 2.790934727908700e-07
% 512 6.977182215317512e-08
