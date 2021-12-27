%clc;
tic;
format short e
clear all;
n=2^9;          
eps=2^(-10);
x(1)=-1;
x(n+1)=1;
h=(x(n+1)-x(1))/n;
y(1)=-1;
y(n+1)=1;
for i=2:n+1
    x(i)=x(1)+(i-1)*h;
end
alf=5/12;   
A2=1/12;
A3=1/12;
omg=-1/(20*eps);
ro=1;
for i=1:n+1
   pp(i)=-2*x(i);
   qq(i)=0;
   ff(i)=0;
end

W(1)=0;
G(1)=-1.000;
for i=2:n 
    
    e(i)=-eps*ro+((A2*h*pp(i+1))/2)-(((alf)*h*pp(i))*(1+(2*omg*h*h*qq(i-1))-(omg*h*(pp(i+1)+3*pp(i-1)))))+(A3*h*h*qq(i-1))-(3*h*A3*pp(i-1)/2);
      f(i)=-(2*eps*ro)+(2*h*A2*pp(i+1))+(4*omg*pp(i)*h*h*(alf)*(pp(i+1)+pp(i-1)))-(qq(i)*h*h*(2*alf))-(2*h*A3*pp(i-1));  
      g(i)=-eps*ro+((3*A2*h*pp(i+1))/2)+(((alf)*h*pp(i))*(1+(2*omg*h*h*qq(i+1))+(omg*h*(3*pp(i+1)+pp(i-1)))))+(A2*h*h*qq(i+1))-((h*A3*pp(i-1))/2);
      hh(i)=-(h*h)*(((A2+(pp(i)*(2*alf)*omg*h))*ff(i+1))+((2*alf)*ff(i))+((A3-((2*alf)*pp(i)*omg*h))*ff(i-1)));  
end

for i=2:n
    W(i)=g(i)/(f(i)-e(i)*W(i-1));
    G(i)=(-hh(i)+e(i)*G(i-1))/(f(i)-e(i)*W(i-1));
end
for i=n:-1:2
    y(i)=G(i)+W(i)*y(i+1);
   end
%y=y';
for i=1:n+1
    y1(i)=erf(x(i)/sqrt(eps));
end
for i=1:n+1
    y2(i)=abs(y(i)-y1(i));
end
y2=y2';

error=abs(y-y1);
max=max(error)'
% format short e
T = [x' y' y1' error'];
%plot(x,y,'+',x,y,'-');
toc

N=[2^(7),2^(8),2^(9),2^(10)];
E1=[7.3536e-06, 4.5818e-07, 2.8613e-08, 1.7890e-09];
E2=[2.9414e-05, 1.8231e-06, 1.1455e-07, 7.1558e-09];
E3=[1.1881e-04, 7.3536e-06, 4.5818e-07, 2.8613e-08];
E4=[4.8395e-04, 2.9414e-05, 1.8231e-06, 1.1455e-07];
%plot(x, y2)
loglog(N, E4,'-*')
legend('\epsilon=2^-^7','\epsilon=2^-^8','\epsilon=2^-^9','\epsilon=2^-^1^0');
xlabel('N-Number of mesh points');
ylabel('Maximum absolute error');  
hold on;

