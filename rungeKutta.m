function [t,sol] = rungeKutta(f,limits,initial,step)

a=[0 0 0 0;
    0.4 0 0 0;
    0.29697761 0.15875964 0 0;
    0.21810040 -3.05096516 3.83286476 0];
b=[0.17476028 -0.55148066 1.20553560 0.17118478];
c=[0 0.4 0.45573725 1];

% a=[0 0 0 0 0 0 0;
%     1/5 0 0 0 0 0 0;
%     3/40 9/40 0 0 0 0 0;
%     44/45 -56/15 32/9 0 0 0 0;
%     19372/6561 -25360/2187 64448/6561 -212/729 0 0 0;
%     9017/3168 -355/33 46732/5247 49/176 -5103/18656 0 0;
%     35/384 0 500/1113 125/192 -2187/6784 11/84 0];
% b=[35/384 0 500/1113 125/192 -2187/6784 11/84 0];
% c=[0 1/5 3/10 4/5 8/9 1 1];

num=floor((limits(end)-limits(1)+1)/step);
sol=zeros(num,length(initial));
sol(1,:)=initial;
for i=1:num-1
    k=zeros(4,length(initial));
    for j=1:4
        k(j,:)=f(limits(1)+step*i+c(j)*step,y(i,:)+step*a(j,:)*k);
    end
    
    if any(isinf(k))
        break
    end
    
    if any(isnan(k))
        break
    end
    
    est=step*b*k;%4th order estimate
    sol(i+1,:)=sol(i,:)+est;%next step of iteration
end
t=step*(1:num);
end