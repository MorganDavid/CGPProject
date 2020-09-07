% Plot basic CGP code
%dataset=table2array(dataset);
% dataset(:,5)=mean(dataset(:,2:4)');
% plot(dataset(:,1),dataset(:,5));
% xlabel("Generation");
% ylabel("Fitness (MSE)");
% title("Simple CGP on 200PointsTrig");
% ax = gca;
% ax.XRuler.Exponent = 0;
% saveas(gca,"basic_cgp200_points_trig.png");

function generate_results_chart()
    harmonicfitness=csvread("harmonicfitness_200pointstrig_waveproperties.csv");
    realfitness=csvread("realfitness_200pointstrig_waveproperties.csv");
    % Remove first row
    harmonicfitness(1,:)=[];
    realfitness(1,:)=[];
    % C++ code doesn't handle early stopped code so fix it here
    harmonicfitness( 1, 2:length(harmonicfitness) - 1 )=0;

    %Remove nans at the end
    realfitness=realfitness(:,1:length(realfitness)-1);
    harmonicfitness=harmonicfitness(:,1:length(harmonicfitness)-1);

    %Flatten the functions
    realfitness=realfitness.';
    realfitness=realfitness(:);
    harmonicfitness=harmonicfitness.';
    harmonicfitness=harmonicfitness(:);
    

    %Plot data when all ready
    x=harmonicfitness.';
    x=x(:);%Change to one column
  % f=figure('visible','on');    
    y=realfitness.';
    y=y(:);%Change to one column
    p=plot(1:20:30000,x,'-',1:20:30000,y,'-','LineWidth',1.2);
    xlabel("Generation");
    ylabel("Fitness (MSE)");
    legend(p([1 2]),"Harmonic MSE","Actual MSE");
    ax = gca;
    ax.XRuler.Exponent = 0;
    xline(00e3,':',{'Epoch 1'},'HandleVisibility','off');
    xline(10e3,':',{'Epoch 2'},'HandleVisibility','off');
    xline(20e3,':',{'Epoch 3'},'HandleVisibility','off');
    title("Harmonic CGP (Approach 2 with Input Approach 2)");
    saveas(gca,"a1_ia2_200pointstrig.png");
end