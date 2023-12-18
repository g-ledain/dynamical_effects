[force, sr] = audioread("guitar.m4a");
force=1.2*force(:,1);
F=griddedInterpolant(1:length(force),force);

freqs=4*[41.2,51.9,61.7,73.4,92.4,110,138.6];%in Hz. This one is Dorian 
freqs=2*pi*freqs;
freqs=freqs/sr;

initial=[0,0];
limits=[1 length(force)];
f= @(t,x,freq) [x(2);-freq^2*B(x(1),1/freq*x(2))-freq^2*A(x(1))+F(t)];
step=1/3;

a=[0 0 0 0;
   0.4 0 0 0;
   0.29697761 0.15875964 0 0;
   0.21810040 -3.05096516 3.83286476 0];
b=[0.17476028 -0.55148066 1.20553560 0.17118478];
c=[0 0.4 0.45573725 1];

sol=zeros((1/step)*length(force),length(initial));
sol(1,:)=initial;

tempsol=zeros((1/step)*length(force),length(initial));
tempsol(1,:)=initial;


for freq=freqs
    tempsol=zeros((1/step)*length(force),length(initial));
    tempsol(1,:)=initial;
    
    for i=1:(1/step)*length(force)-1
        k=zeros(4,length(initial));
        for j=1:4
            k(j,:)=f(limits(1)+step*i+c(j)*step,tempsol(i,:)+step*a(j,:)*k,freq);
        end

        if any(any(isinf(k))) || any(any(isnan(k)))
            break
        end

        tempsol(i+1,:)=tempsol(i,:)+step*b*k;%next step of iteration
    end
    sol=sol+tempsol;
end
sol=sol./length(freqs);
x=sol(:,1);
x=x(1:(1/step):end);
audiowrite("manyOsc.wav",x/max(abs(x)),sr)
%combo=([force;0]+x)/max(abs([force;0]+x));
%audiowrite("force.wav",force,sr)
%audiowrite("combomanyosc.wav",([force;0]+x)/max(abs([force;0]+x)),sr)


%%
function y=A(x)
y=x;
end

function y=B(x,v)
y=-(1-x^2)*v;
end