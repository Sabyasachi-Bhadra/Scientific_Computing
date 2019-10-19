a=0;b=1;
%dividing this interval in N equal parts
N=16;
%h be the steplength
h=(b-a)/N;
x=linspace(a,b,N+1);

%taking the transpose of x
x=x';

%r.h.s of the given differential equation
f=exp(2*x);

%excluding the two boundary condition
X=x(2:N);
%same for the r.h.s
F=exp(2*X);

%produce the sparse matrix A
n=N-1;e=ones(n,1);
A=((-1)/(h*h))*spdiags([e -2*e e],-1:1,n,n);
sparse(A)
%full(A)



%solving the equation AU=F
U=A\F;



%analytical solution of u
u=((-1)/4)*(exp(2*X)+(1-exp(2))*X-1);
u=[0;u;0];
U=[0;U;0];




%difference between U and u




diff=max(abs(U-u));
%plotting in the diagram
plot(x,U,'-o')
hold on
plot(x,u,'-*')
legend('computational','analytical')
