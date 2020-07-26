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
% Abar = [A zeros(6,2);-C eye(2,2)];
% Bbar = [B;zeros(2,2)];
% rank(ctrb(Abar,Bbar))
desired_poles = [-1+1j -1-1j -2+2j -2-2j -2 -2];
K = place(A,B,desired_poles);
G0 = -C/(A-B*K)*B;
%% LINEAR ODE45
init = [0 0 5*3.14/180 -5*3.14/180 10*3.14/180 -10*3.14/180];
% init = [0 0 5*3.14/180 5*3.14/180 2*3.14/180 -1*3.14/180];
options = odeset('RelTol',1e-2,'AbsTol',1e-4);
tspan = 0:0.01:100;
[t,X] = ode45(@(t,x) linear_ode(t,x,A,B,G0,K),tspan,init,options);
subplot(3,2,1);plot(t,X(:,1));title('X');
subplot(3,2,2);plot(t,X(:,2));title('X-dot');
subplot(3,2,3);plot(t,X(:,3)/3.14*180);title('Teta1');
subplot(3,2,4);plot(t,X(:,4)/3.14*180);title('Teta-dot1');
subplot(3,2,5);plot(t,X(:,5)/3.14*180);title('Teta2');
subplot(3,2,6);plot(t,X(:,6)/3.14*180);title('Teta-dot2');
figure;
[t,XX] = ode45(@(t,x) nonlinear_ode(t,x,G0,K),tspan,init,options);
subplot(3,2,1);plot(t,XX(:,1));title('X');
subplot(3,2,2);plot(t,XX(:,2));title('X-dot');
subplot(3,2,3);plot(t,XX(:,3)/3.14*180);title('Teta1');
subplot(3,2,4);plot(t,XX(:,4)/3.14*180);title('Teta-dot1');
subplot(3,2,5);plot(t,XX(:,5)/3.14*180);title('Teta2');
subplot(3,2,6);plot(t,XX(:,6)/3.14*180);title('Teta-dot2');

Yr1 = 0.3*sign(sin(0.2*t));
Yr2 = 0.3*sign(sin(0.2*t))*3.1415/180;
figure
plot(t,Yr1,t,X(:,1))
figure
plot(t,Yr2,t,X(:,3))
figure
plot(t,Yr1,t,XX(:,1))
figure
plot(t,Yr2,t,XX(:,3))



function dx = nonlinear_ode(t,x,G0,K)
    q1=x(1);dq1=x(2);
    q2=x(3);dq2=x(4);
    q3=x(5);dq3=x(6);

    Yr = 0.3*sign(sin(0.2*t))*[1;3.1415/180];
    u = -K*x+G0\Yr;u1 = u(1);u2 = u(2);
   
    ddq =[(1962*sin(2*q2) - 1500*u1 - 200*dq2^2*sin(q2) + 100*dq2^2*sin(q3) + 500*u1*cos(2*q2 - 2*q3) + 50*dq2^2*cos(2*q2 - 2*q3) - 75*dq2^2*sin(2*q2 - q3) + 25*dq2^2*sin(2*q2 - 3*q3) + 3000*u2*cos(q2) - 150*dq2^2 - 1000*u2*cos(q2 - 2*q3))/(200*(cos(2*q2) + 5*cos(2*q2 - 2*q3) - 16));
          -(23000*u2 + 9810*sin(q2 - 2*q3) + 33354*sin(q2) - 1000*u2*cos(2*q3) + 50*dq2^2*cos(q2 - 2*q3) - 1150*dq2^2*sin(q2 - q3) + 25*dq2^2*sin(q2 - 3*q3) - 100*dq2^2*sin(2*q2) - 500*dq2^2*sin(2*q2 - 2*q3) - 1500*u1*cos(q2) + 25*dq2^2*sin(q2 + q3) + 500*u1*cos(q2 - 2*q3) - 150*dq2^2*cos(q2))/(100*(cos(2*q2) + 5*cos(2*q2 - 2*q3) - 16));
          -(1962*sin(q3) - 1962*sin(2*q2 - q3) + 200*dq2^2*sin(q2 - q3) + 100*u1*cos(2*q2 - q3) - 5*dq2^2*sin(2*q2) + 5*dq2^2*sin(2*q3) + 200*u2*cos(q2 + q3) + 10*dq2^2*cos(2*q2 - q3) + 55*dq2^2*sin(2*q2 - 2*q3) - 100*u1*cos(q3) - 2200*u2*cos(q2 - q3) - 10*dq2^2*cos(q3))/(10*(cos(2*q2) + 5*cos(2*q2 - 2*q3) - 16))];

    dx(1) = dq1;
    dx(2) = ddq(1);
    dx(3) = dq2;
    dx(4) = ddq(2);
    dx(5) = dq3;
    dx(6) = ddq(3);
    dx=[dx(1);dx(2);dx(3);dx(4);dx(5);dx(6)];
end
function dX = linear_ode(t,X,A,B,G0,K)
    Yr = 0.3*sign(sin(0.2*t))*[1;3.1415/180];
    u = -K*X+G0\Yr;
    dX = A*X + B*u;
end
