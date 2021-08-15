 %Resolution of Blasius  equation(with Runge-Kutta)
clc
clear all
a=0;
b=10;
h=0.02;
N=(b-a)/h;
F2=zeros(91,N+1);
k=0.3321; %f''(0)
y=[0;0;k];
y1=y(1);%f
y2=y(2);%f'
y3=y(3);%f''
for n=1:N
    A=[0 1 0;0 0 1;(-y(3)/2) 0 0];
    k1=A*y;
    k2=A*(y+(h/2)*k1);
    k3=A*(y+(h/2)*k2);
    k4=A*(y+h*k3);
    y=y+(h/6)*(k1+2*k2+2*k3+k4);
    y1(n+1)=y(1);
    y2(n+1)=y(2);
    y3(n+1)=y(3);
end
t=[a:h:b];
num=int8(100*(k-0.1+0.01));
%F2 has as rows approximations of y2 for the different values of k
F2(num,:)=y2;
grid on
t=[a:h:b];
k1=[0.1:0.01:1];
f10=F2(:,N+1);
f10=f10';
o=ones(1,91);
bigmatrix=[t.',y1.',y2.',y3.']%table of eta f f' f''
figure(1)
hold on
plot(y2,t,'b')
ylabel('\eta')
xlabel( 'f´ or u/U')
hold off