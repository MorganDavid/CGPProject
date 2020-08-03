csvfile = csvread("plot_1.csv");
f = figure('visible','off');
plot(csvfile(:,1),csvfile(:,2),csvfile(:,1),csvfile(:,3));
legend("true","pred");
saveas(f,"plot_1.jpg");

csvfile = csvread("plot_2.csv");
f = figure('visible','off');
plot(csvfile(:,1),csvfile(:,2),csvfile(:,1),csvfile(:,3));
legend("pred","true");
saveas(f,"plot_2.jpg");

csvfile = csvread("plot_final.csv");
f = figure('visible','off');
plot(csvfile(:,1),csvfile(:,2),csvfile(:,1),csvfile(:,3));
legend("pred","true");
saveas(f,"plot_final.jpg");

csvfile = csvread("harmonics.csv.csv");
sz = size(csvfile);
f = figure('visible','off');
hold on 
for x = 1:sz(1)
    plot(csvfile(x,:),"DisplayName",num2str(x)+" harm(s)");
end
hold off
legend show
saveas(f,"harmonics.jpg");
