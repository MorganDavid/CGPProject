%To plot the actual vs harmonic fitness diagram. 

t=1:0.01:4;
x_actual =  4*sin(2*pi*4*t+1)+sin(2*pi*3*t+0.3)+cos(2*pi*9*t+0.2)*2;
x_prediction = 3.8*sin(2*pi*4*t+1.1);
x_harmonic =  4*sin(2*pi*4*t+1);
plot(t,x_actual,t,x_harmonic,t,x_prediction,'LineWidth',1.2);
legend("True Function ($Y$)","Synthesis at epoch ($Z_1$)","Prediction at epoch 1 ($\hat{Y}_{1}$)",'Interpreter','latex','FontSize',12)
%actual MSE: 2.5983, harmonic MSE: 0.0958
set(gca,'FontSize',12)
saveas(gcf,'actual_vs_harmonic.png')
