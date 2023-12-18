[force, sr] = audioread("guitar.m4a");
force=force(82569:end,1);%trim audio
force=10.*force;
F=griddedInterpolant(1:length(force),force);

a=[0 0 0 0;
   0.4 0 0 0;
   0.29697761 0.15875964 0 0;
   0.21810040 -3.05096516 3.83286476 0];
b=[0.17476028 -0.55148066 1.20553560 0.17118478];
c=[0 0.4 0.45573725 1];


initial=[0,0];
limits=[1 length(force)];

K=1/(3*pi)*22050/20000;
f= @(t,x) [x(2);K*9*(1-x(1)^2)*x(2)-K^2*x(1)+K^2*F(t)];
%f= @(t,x) [x(2);K*9*(1-x(1)^2)*x(2)+K^2*F(t)];
step=1;

[t,sol] = rungeKutta(f,limits,initial,step)

audiowrite("forcedvdp.wav",sol(:,1)/max(abs(sol(:,1))),sr)
%audiowrite("force.wav",force,sr)
audiowrite("combovdp.wav",([force;0]+sol(:,1))/max(abs([force;0]+sol(:,1))),sr)
