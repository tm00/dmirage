v=0.005;
tf=200;
time_steps=200;
delta=tf/time_steps;
bf=0.5;

% Bfield1
for n=0:tf
    t = n*delta
    energy(n+1)=sum(Eprof(:,n+1)) + (1-v*t)*bf*S1prof(10,n+1)+v*t*bf*S1prof(11,n+1);
end
