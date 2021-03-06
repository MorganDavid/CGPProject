dataset=basiccgp200pointstrigresults;
%dataset=table2array(dataset);
dataset(:,5)=mean(dataset(:,2:4)');
plot(dataset(:,1),dataset(:,5));
xlabel("Generation");
ylabel("Fitness (MSE)");
title("Simple CGP on 200PointsTrig");
ax = gca;
ax.XRuler.Exponent = 0;
saveas(gca,"basic_cgp200_points_trig.png");