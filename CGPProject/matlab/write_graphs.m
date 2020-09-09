% Writes all the results predictions and harmonics outputted into the /output folder. 

subdir="images/";

csvfile = csvread("plot_1.csv");
f = figure('visible','off','PaperPosition',[0 0 38 11]);
title("One harmonics pred vs. true");
plot(csvfile(:,1),csvfile(:,2),csvfile(:,1),csvfile(:,3));
legend("pred","true");
saveas(f,subdir+"plot_1.jpg");

csvfile = csvread("plot_2.csv");
f = figure('visible','off','PaperPosition',[0 0 38 11]);
title("Two harmonics pred vs. true");
plot(csvfile(:,1),csvfile(:,2),csvfile(:,1),csvfile(:,3));
legend("pred","true");
saveas(f,subdir+"plot_2.jpg");

csvfile = csvread("plot_final.csv");
f = figure('visible','off','PaperPosition',[0 0 38 11]);
title("Three harmonics pred vs. true");
plot(csvfile(:,1),csvfile(:,2),csvfile(:,1),csvfile(:,3));
legend("pred","true");
saveas(f,subdir+"plot_final.jpg");

csvfile = csvread("harmonics.csv");
sz = size(csvfile);
f = figure('visible','off','PaperPosition',[0 0 80 20]);
title("Synthesised harmonics from Fourier transform");
hold on 
for x = 1:sz(1)
    plot(csvfile(x,:),"DisplayName",num2str(x)+" harm(s)");
end
hold off
legend show
saveas(f,subdir+"harmonics.tiff");

% z = peaks(25);
% figure
% mesh(z)
% xlabel("x");
% ylabel("y");
% zlabel("z");
% text(13,19,8+0.4,'Global Optimum',"FontSize",12);
% text(18,13,3.5+0.6,'Local Optimum',"FontSize",12);
% set(gca,'FontSize',12)
