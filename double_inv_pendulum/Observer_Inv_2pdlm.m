clc;
clear;
close all;

global M m1 m2 L1 L2 g d1 d2 d3 w1 w2 w3 
M = 2;m1 = 0.2;m2 = m1;
L1 = 0.5;L2 = L1;
g = 9.81;
d1=0;d2=0;d3=0;
w1=0;w2=0;w3=0;

[A,B,C,D]=State_Space();
%%
rank(ctrb(A,B))
desired_poles = [-1+1j -1-1j -2+2j -2-2j -2 -2];
K = place(A,B,desired_poles);

%% EESA METHOD
% landa = eig(A);
% V = zeros(6,6);q = zeros(2,6);
% for i=1:6
%     OB = [A'-landa(i)*eye(6,6) C'];
%     N = null(OB);
%     V(:,i) = N(1:6,randi([1 2],1,1));
%     q(:,i) = N(7:8,randi([1 2],1,1));
% end
% % L = (-q/V)';
% L=[9.0236   59.2801
%    -0.0000  367.5931
%    -0.8484  -21.4125
%     8.7613  160.5018
%     0.5498   24.4832
%   -17.5226 -235.5249];
%% POLE PLACEMENT METHOD
rank(obsv(A,C))
des_poles = [-5.01 -5 -10.01 -10 -15.01 -15];
L = place(A',C',des_poles)';
%% LINEAR ODE45
% init = 0.1*[0 0 5*3.14/180 -5*3.14/180 10*3.14/180 -10*3.14/180 zeros(1,6)];
init = 0.1*[0 0 5*3.14/180 5*3.14/180 2*3.14/180 -1*3.14/180 zeros(1,6)];
% init = [0 0 1*3.14/180 0.1*3.14/180 1*3.14/180 -1*3.14/180 zeros(1,6)];
options = odeset('RelTol',1e-2,'AbsTol',1e-4);
tspan = 0:0.01:10;
% [t,X] = ode45(@(t,x) linear_ode(t,x,A,B,C,K,L),tspan,init,options);
% subplot(3,2,1);plot(t,X(:,1));title('X');
% subplot(3,2,2);plot(t,X(:,2));title('X-dot');
% subplot(3,2,3);plot(t,X(:,3)/3.14*180);title('Teta1');
% subplot(3,2,4);plot(t,X(:,4)/3.14*180);title('Teta-dot1');
% subplot(3,2,5);plot(t,X(:,5)/3.14*180);title('Teta2');
% subplot(3,2,6);plot(t,X(:,6)/3.14*180);title('Teta-dot2');
% figure;
[t,XX] = ode45(@(t,x) nonlinear_ode(t,x,K,L,C),tspan,init,options);
subplot(3,2,1);plot(t,XX(:,1));title('X');
subplot(3,2,2);plot(t,XX(:,2));title('X-dot');
subplot(3,2,3);plot(t,XX(:,3)/3.14*180);title('Teta1');
subplot(3,2,4);plot(t,XX(:,4)/3.14*180);title('Teta-dot1');
subplot(3,2,5);plot(t,XX(:,5)/3.14*180);title('Teta2');
subplot(3,2,6);plot(t,XX(:,6)/3.14*180);title('Teta-dot2');
figure;
%% PLOT
% plot(t,X(:,1),t,XX(:,1),'r')
% legend('Linear Equation','NonLinear Equation')
% title('Compare X')
% xlabel('Time');ylabel('X')
% figure;
% plot(t,X(:,3),t,XX(:,3),'r')
% legend('Linear Equation','NonLinear Equation')
% title('Compare Teta')
% xlabel('Time');ylabel('Teta')
% figure
% plot(t,X(:,1),t,X(:,7),'r')
% legend('without OBSERVER','with OBSERVER')
% title('Compare X')
% xlabel('Time');ylabel('X')
% figure;
% plot(t,X(:,3),t,X(:,9),'r')
% legend('without OBSERVER','with OBSERVER')
% title('Compare Teta1')
% xlabel('Time');ylabel('Teta1')
% figure;
% plot(t,X(:,5),t,X(:,11),'r')
% legend('without OBSERVER','with OBSERVER')
% title('Compare Teta2')
% xlabel('Time');ylabel('Teta2')
% figure;
plot(t,XX(:,1),t,XX(:,7),'r')
legend('without OBSERVER','with OBSERVER')
title('Compare X')
xlabel('Time');ylabel('X')
figure;
plot(t,XX(:,3),t,XX(:,9),'r')
legend('without OBSERVER','with OBSERVER')
title('Compare Teta1')
xlabel('Time');ylabel('Teta1')
figure;
plot(t,XX(:,5),t,XX(:,11),'r')
legend('without OBSERVER','with OBSERVER')
title('Compare Teta2')
xlabel('Time');ylabel('Teta2')
%%
function dX = nonlinear_ode(t,x,K,L,C)
    q1=x(1);dq1=x(2);
    q2=x(3);dq2=x(4);
    q3=x(5);dq3=x(6);
    u = -K*x(7:12);u1 = u(1);u2 = u(2);
    
    Y = C*x(1:6);
    Yh = C*x(7:12);
   
    ddq =[u1/2 - u2 - (981*sin(q2))/500 + (dq2^2*sin(q2))/10 + (dq2^2*sin(q2 - q3))/20 + dq2^2/20;
          22*u2 - u1 + (10791*sin(q2))/250 - (981*sin(q3))/50 - (dq2^2*sin(q2))/5 - (21*dq2^2*sin(q2 - q3))/10 - dq2^2/10;
          3*sin(q2 - q3)*dq2^2 - 20*u2 - (981*sin(q2))/25 + (981*sin(q3))/25];
    dX(1) = dq1;
    dX(2) = ddq(1);
    dX(3) = dq2;
    dX(4) = ddq(2);
    dX(5) = dq3;
    dX(6) = ddq(3);
    dX = [dX(1);dX(2);dX(3);dX(4);dX(5);dX(6)];
    
    qh1=x(7);dqh1=x(8);
    qh2=x(9);dqh2=x(10);
    qh3=x(11);dqh3=x(12);
    ddqh =[u1/2 - u2 - (981*sin(qh2))/500 + (dqh2^2*sin(qh2))/10 + (dqh2^2*sin(qh2 - qh3))/20 + dqh2^2/20;
          22*u2 - u1 + (10791*sin(qh2))/250 - (981*sin(qh3))/50 - (dqh2^2*sin(qh2))/5 - (21*dqh2^2*sin(qh2 - qh3))/10 - dqh2^2/10;
          3*sin(qh2 - qh3)*dqh2^2 - 20*u2 - (981*sin(qh2))/25 + (981*sin(qh3))/25];
    dXh(1) = dqh1;
    dXh(2) = ddqh(1);
    dXh(3) = dqh2;
    dXh(4) = ddqh(2);
    dXh(5) = dqh3;
    dXh(6) = ddqh(3);
    dXh = [dXh(1);dXh(2);dXh(3);dXh(4);dXh(5);dXh(6)];
    dXh = dXh + L*(Y-Yh); 
    dX=[dX;dXh];
end
function dx = linear_ode(t,XX,A,B,C,K,L)
    X = XX(1:6);
    Xhat=XX(7:12);
    u = -K*Xhat;
    dX = A*X + B*u;
    Y = C*X;
    Yhat = C*Xhat;
    dXhat = A*Xhat + B*u + L*(Y-Yhat);
    dx=[dX;dXhat];
end
