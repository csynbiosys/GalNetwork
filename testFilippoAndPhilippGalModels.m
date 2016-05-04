function testFilippoAndPhilippGalModels
% This code refers to the committ: 543052642858a811a8f4ec394b36e479f47c145e

[t1,y1]=ode45(@montiGALphilipp,[0 3000],[31.9700    0.0001    0.3630   34.1947]);
[t2,y2]=ode45(@montiGALmodel,[0 3000],[31.9700    0.0001    0.3630   34.1947]);

plot(t1,y1)
hold on
plot(t2,y2)
hold off
