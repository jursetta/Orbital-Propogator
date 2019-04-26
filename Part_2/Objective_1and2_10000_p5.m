dataRaw1 = dlmread('Optimum_1_10000_p5');
dataRaw2 = dlmread('Optimum_2_10000_p5');

r_E = 6371000;
ex  = linspace(-r_E,r_E);
eyu =  sqrt(-ex.^2+r_E^2);
eyl = -sqrt(-ex.^2+r_E^2);

figure
title('Flight Path of objective 1 for 10,000(m) clearance and .5(m/s) \Delta V Tol')
xlabel('Distance (m)')
ylabel('Distance (m)')
hold on
plot(dataRaw1(:,2) ,dataRaw1(:,3))
plot(dataRaw1(:,6) ,dataRaw1(:,7))
plot(dataRaw1(:,10),dataRaw1(:,11))
legend('Space Craft Trajectory','Lunar Trajectory','Location','SouthEast')
plot(ex,eyu,'k')
plot(ex,eyl,'k')
axis equal
saveas(gcf,'Objective1.jpg')

figure
title('Flight Path of objective 2 for 10,000(m) clearance and .5(m/s) \Delta V Tol')
xlabel('Distance (m)')
ylabel('Distance (m)')
hold on
plot(dataRaw2(:,2) ,dataRaw2(:,3))
plot(dataRaw2(:,6) ,dataRaw2(:,7))
plot(dataRaw2(:,10),dataRaw2(:,11))
legend('Space Craft Trajectory','Lunar Trajectory','Location','SouthEast')
plot(ex,eyu,'k')
plot(ex,eyl,'k')
axis equal
saveas(gcf,'Objective2.jpg')