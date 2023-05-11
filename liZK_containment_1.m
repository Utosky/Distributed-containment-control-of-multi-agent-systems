%paper title:   Distributed containment control of multi-agent systems with general linear dynamics in the presence of multiple leaders
%paper author:  Zhongkui Li 2011
%example 1:     continuous-time static containment controller
%code author:   Zhiwen.He
 
%% initial state
clc;clear;close all;
global A B L1 L2
A=[0 0    0       1      0      0;
   0 0    0       0      1      0;
   0 0    0       0      0      1;
   0 0 -0.2003 -0.2003   0      0;
   0 0  0.2003    0   -0.2003   0;
   0 0    0       0      0    -1.6129];
B=[ 0       0;
    0       0;
    0       0;
  0.9441 0.9441;
  0.9441 0.9441;
-28.7079 28.7079];

L1=[ 3  0  0 -1 -1 -1 ;
   -1  1  0  0  0  0 ;
   -1 -1  2  0  0  0 ;
   -1  0  0  2  0  0 ;
    0  0  0 -1  2  0;
    0  0  0  0 -1  2];
L2=[0  0  0;
    0  0  0;
    0  0  0;
    0  0 -1;
    0 -1  0;
   -1  0  0];

L=[ 3  0  0 -1 -1 -1  0  0  0;
   -1  1  0  0  0  0  0  0  0;
   -1 -1  2  0  0  0  0  0  0;
   -1  0  0  2  0  0  0  0 -1;
    0  0  0 -1  2  0  0 -1  0;
    0  0  0  0 -1  2 -1  0  0;
    0  0  0  0  0  0  0  0  0;
    0  0  0  0  0  0  0  0  0;
    0  0  0  0  0  0  0  0  0];
[V,D]=eig(L1);
c_th=1/min(diag(D));        % c_th=1.2176, c>=c_th,取c=2

%% LMI求解P
setlmis([])
P=lmivar(1,[6 1]);
lmiterm([1 1 1 P],1,A','s');
lmiterm([-1 1 1 0],2*B*B');
lmiterm([-2 1 1 P],1,1);
lmis1 = getlmis;
[tmin,xfeas] = feasp(lmis1);
%判断是否有可行解
if (tmin < 0)  
    disp('Feasible');
    P = dec2mat(lmis1,xfeas,P);
else
    P = nan;
end
global F
F=-B'*inv(P);

% F=[-0.0089, 0.0068, 0.0389,-0.0329, 0.0180, 0.0538;
%     0.0068,-0.0089,-0.0389, 0.0180,-0.0329,-0.0538];

%% calculate ODE funtion
tspan=[0,50];
x=2*rand(54,1);
[t,y]=ode45(@contain,tspan,x);

%% draw picture

figure(1)
% plot(t,y(:,1),t,y(:,7),t,y(:,13),t,y(:,19),t,y(:,25),t,y(:,31),t,y(:,37),t,y(:,43),t,y(:,49),'linewidth',1.5);
for i=1:8
    if(i<6)
        plot(t,y(:,1+6*i),'b','linestyle',':','linewidth',1.5);
    else
        plot(t,y(:,1+6*i),'r','linewidth',2);
    end
grid on;
hold on;
end
xlabel('t');
ylabel('x_{i1}');
% figure(1)
% plot(t,y(:,1),t,y(:,7),t,y(:,13),t,y(:,19),t,y(:,25),t,y(:,31),'linewidth','0.8');
% % plot(t,y(:,37),t,y(:,43),t,y(:,49),'linewidth',3.5);
% plpt(t,y(:,37),t,y(:,43),t,y(:,49),'linewidth','3.5','color','r');
xlabel('t');
ylabel('x_{i1}');

figure(2)
for i=1:8
    if(i<6)
        plot(t,y(:,2+6*i),'b','linestyle',':','linewidth',1);
    else
        plot(t,y(:,2+6*i),'r','linewidth',1.5);
    end
grid on;
hold on;
end
xlabel('t');
ylabel('x_{i2}');

% figure(2)
% for i=1:8
% plot(t,y(:,2+6*i),'LineWidth',1.5);
% grid on;
% hold on;
% end
% xlabel('t');
% ylabel('x_{i2}');

% figure(3)
% for i=1:8
% plot(t,y(:,3+6*i),'LineWidth',1.5);
% grid on;
% hold on;
% end
% legend('$x3_1$','$x3_2$','$x3_3$','$x3_4$','$x3_5$','$x3_6$','$x3_7$','$x3_8$','$x3_9$','Interpreter','latex');
% xlabel('t');
% ylabel('x_{i3}');

figure(3)
for i=1:8
    if(i<6)
        plot(t,y(:,3+6*i),'b','linestyle',':','linewidth',1);
    else
        plot(t,y(:,3+6*i),'r','linewidth',1.5);
    end
grid on;
hold on;
end
xlabel('t');
ylabel('x_{i3}');


%% function 
function   ydot=contain(t,x)
c=2;  
global A B F L1 L2
a11=kron(eye(6),A)+c*kron(L1,B*F);
a12=c*kron(L2,B*F);
a22=kron(eye(3),A);
a21=zeros(18,36);
ydot=[a11 a12;
      a21 a22]*x;
end

