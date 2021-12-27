tic%clc;
clear all;
n=2^10;          
eps=2^(-10);
%eta=0.1*eps;
x(1)=0;
x(n+1)=1;
h=(x(n+1)-x(1))/n;
y(1)=1;
y(n+1)=1;
for i=2:n+1
    x(i)=x(1)+(i-1)*h;
end
alf=5/12;
A2=1/12;
A3=1/12;
omg=-1/(20*eps);

%d=0.1*eps;

ro=1;

for i=1:n+1
   pp(i)=-2*(2*x(i)-1);
   qq(i)=4;%1)2000./(t+1);
   ff(i)=0;%(exp(x(i)*x(i))*(2*eps+4*eps*x(i)*x(i)-3-x(i)*x(i)));
%(2*eps-1)*cosh(x(i))+eps*x(i)*sinh(x(i))-sinh(x(i))/x(i);
end

W(1)=0;
G(1)=1.000;
% e(2)=-1+((alf*h*h*q(1))/(eps-k));
% f(2)=-2-((2*bet*(h*h)*q(2))/(eps-k));
% g(2)=-1+((alf*h*h*q(3))/(eps-k));
% hh(2)=-((h*h)*(alf*c(1)+2*bet*c(2)+alf*c(3)))/(eps-k);
for i=2:n
    e(i)=-eps*ro+((A2*h*pp(i+1))/2)-(((alf)*h*pp(i))*(1+(2*omg*h*h*qq(i-1))-(omg*h*(pp(i+1)+3*pp(i-1)))))+(A3*h*h*qq(i-1))-(3*h*A3*pp(i-1)/2);
      f(i)=-(2*eps*ro)+(2*h*A2*pp(i+1))+(4*omg*pp(i)*h*h*(alf)*(pp(i+1)+pp(i-1)))-(qq(i)*h*h*(2*alf))-(2*h*A3*pp(i-1));  
      g(i)=-eps*ro+((3*A2*h*pp(i+1))/2)+(((alf)*h*pp(i))*(1+(2*omg*h*h*qq(i+1))+(omg*h*(3*pp(i+1)+pp(i-1)))))+(A2*h*h*qq(i+1))-((h*A3*pp(i-1))/2);
      hh(i)=-(h*h)*(((A2+(pp(i)*(2*alf)*omg*h))*ff(i+1))+((2*alf)*ff(i))+((A3-((2*alf)*pp(i)*omg*h))*ff(i-1)));  
%(1-h*(1-0.001/h)*p(i+1)/2-h*q(i+1)/2);
% 
%       e(i)=-eps-(3*alf*h*p(i-1))/2+(alf*(h*h)*q(i-1))-(h*bet*p(i))-(2*(h^3)*bet*omg*p(i)*q(i-1))+(h*h*bet*omg*p(i)*(p(i+1)+3*p(i-1)))+(alf*h*p(i+1))/2;
%       f(i)=-(2*eps)+(2*alf*h*p(i+1))+(4*h*h*bet*omg*p(i)*(p(i+1)+p(i-1)))-(2*h*h*bet*q(i))-(2*alf*h*p(i-1));  
%       g(i)=-eps+(3*alf*h*p(i+1))/2+(alf*(h*h)*q(i+1))+(h*bet*p(i))+(2*(h^3)*bet*omg*p(i)*q(i+1))+(h*h*bet*omg*p(i)*(3*p(i+1)+p(i-1)))-(alf*h*p(i-1))/2;
%       hh(i)=-(h*h)*((alf+(2*bet*h*omg*p(i)))*c(i+1)+(2*bet)*c(i)+(alf-(2*bet*h*omg*p(i)))*c(i-1));  
% %(1-h*(1-0.001/h)*p(i+1)/2-h*q(i+1)/2);
end

for i=2:n
    W(i)=g(i)/(f(i)-e(i)*W(i-1));
    G(i)=(-hh(i)+e(i)*G(i-1))/(f(i)-e(i)*W(i-1));
end
for i=n:-1:2
    y(i)=G(i)+W(i)*y(i+1);
    %y(n)=1/sqrt(exp(1));
    % y(i+1)=(hh(i)-e(i)*y(i-1)+f(i)*y(i))/g(i);
end
%y=y';
for i=1:n+1
   A(i)=exp(-(2*x(i)*x(i)-2*x(i))/eps);%exp(1/(2*eps)-(1-2*x(i))*(1-2*x(i))/(2*eps));% y(k)=i*(i+1-.002)+(0.002-1)*(1-exp(-1000*i))/(1-exp(-1000));
   B(i)=exp((1-2*x(i))*(1-2*x(i))/(2*eps)); %y(k)=(exp((1.0001*(i-1))*10000))+exp(-i);
   C(i)=sqrt(2*pi)*erf((1-2*x(i))/sqrt(2*eps));       % y1(i)=((exp(c2)-1)*exp(c1*x(i))+(1-exp(c1))*exp(c2*x(i)))/(exp(c2)-exp(c1));
   D(i)=exp(1/(2*eps))*sqrt(2*pi)*erf(1/sqrt(2*eps))+2*sqrt(eps);     %y1(i)=2*exp(-2*(x(i)+1)/eps)+exp(-2*(1-x(i))/eps);
%yy(i)=-A(i)*(2*B(i)*C(i)*x(i)-B(i)*C(i)-2*sqrt(eps))/(D(i));
%yy(i)=-(exp(1/(2*eps)-(1-2*x(i))*(1-2*x(i))/(2*eps)))*(2*(exp((1-2*x(i))*(1-2*x(i))/(2*eps)))*(sqrt(2*pi)*erf((1-2*x(i))/sqrt(2*eps)))*x(i)...
  %  -(exp((1-2*x(i))*(1-2*x(i))/(2*eps)))*(sqrt(2*pi)*erf((1-2*x(i))/sqrt(2*eps)))-2*sqrt(eps))/(exp(1/(2*eps))*sqrt(2*pi)*erf(1/sqrt(2*eps))+2*sqrt(eps));
y1(i)=-(2*C(i)*x(i)-C(i)-2*sqrt(eps)*exp(-(1-2*x(i))*(1-2*x(i))/(2*eps)))/(sqrt(2*pi)*erf(1/sqrt(2*eps))+2*sqrt(eps)*exp(-1/(2*eps)));

end
%    y1(i)=exp(x(i)*x(i));%x(i)*sinh(x(i));
%y1=y1';
format short e
error=abs(y-y1);
max=max(error)'
% format short e
T = [x' y' y1' error'];%zh1=abs((y-y1)/abs(y))';
for i=1:n+1
    y2(i)=abs(y(i)-y1(i));
end
y2=y2';
toc%maxh1=max(zh1)';
N=[2^(7),2^(8),2^(9),2^(10)];
E1=[2.7061e-07, 1.6919e-08, 1.0576e-09, 6.6184e-11];
E2=[7.6654e-07, 4.7874e-08, 2.9916e-09, 1.8699e-10];
E3=[2.1700e-06, 1.3531e-07, 8.4593e-09, 5.2878e-10];
E4=[6.1444e-06, 3.8327e-07, 2.3937e-08, 1.4958e-9];
%plot(x, y2)
loglog(N, E4,'-*')
legend('\epsilon=2^-^7','\epsilon=2^-^8','\epsilon=2^-^9','\epsilon=2^-^1^0');
xlabel('N-Number of mesh points');
ylabel('Maximum absolute error');  
hold on;
%M=max(maxh1)
%plot(x,y1,'-',x,y,'.')
%plot(x,y2,'-');
%hold on;
%legend(num2str('exact solution'),num2str('numerical solution'))% order=(log(maxh1)-log(maxh2))/log(2);
%E1=sqrt(sum(error))/sqrt(sum(y1'*y1))

%c1=(-1-sqrt(1+4*(eps-d)))/(2*(eps-d));
% c2=(-1+sqrt(1+4*(eps-d)))/(2*(eps-d));
%  for i=1:n+1
%  ye(i)=((1-exp(c2))*exp(c1*x(i))+(exp(c1)-1)*exp(c2*x(i)))/(exp(c1)-exp(c2));
% 
% % % %y1(i)=2*x(i)+(1-exp(-x(i)/eps))/(exp(-1/eps)-1);
%  end
% ye=ye';
% maxh2=max(zh2)
% order=(log(maxh1)-log(maxh2))/log(2);
% for i=1:n
%     ff(i)=e(i)+g(i);
% end
% TT=[ff' f'];
% ss=[e1' g1'];
% plot(x,y,'-');
% 
% 
% 
% k=1;
% for i=1:2:n1+2
%     x11(k)=x1(i);
%     y11(k)=y1(i);
%     k=k+1;
%     end
% y11=y11';
% k=1;
% for i=1:4:n2+2
%     x22(k)=x2(i);
%     y22(k)=y2(i);
%     k=k+1;
% end
% y22=y22';
% 
% T=[x' y y11 y22];
% zh1=abs(y-y11);
% zh2=abs(y11-y22);
% maxh1=max(zh1)
% maxh2=max(zh2)
% order=(log(maxh1)-log(maxh2))/log(2);
% for i=1:n
%     ff(i)=e(i)+g(i);
% end
% TT=[ff' f'];
% ss=[e1' g1'];
 %plot(x,ye,'-');
 %plot(x,ye,'-',x,y,'o');
% 
