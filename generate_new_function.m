f=4
t=-1:0.01:1
n=length(t)
%Func1 y=2*sin(2*pi*f*t)+cos(2*pi*1*t)*3;
%Func2 y=sin(t)+sin(t+t.^2);
%Func3 y=t.^2+t.^2+t;
y= t.^4+t.^3+t.^2+t;
plot(t,y)
numberOfSamplesToTake = 50;
sampleIndexes = randperm(numel(y), numberOfSamplesToTake)
% Plot the samples;
ts = t(sampleIndexes)
ys = y(sampleIndexes)
hold on;
plot(ts, ys,'*');

func4=horzcat(ts.',ys.');