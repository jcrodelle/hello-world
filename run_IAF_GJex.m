% This script calculates, for all models, the coupling coeff.
clear all
% close all
% parameters used in all sims:
C = 1.0; 
gL = 0.1; 
Eleak = -70; %mV
vR = -65;
%vE = -90.0;
vE = 0.0;
vT = -50; % for expIAF
sigmar = 1.5;
s = 0.0; % 0.25; %0.09;
%intial conditions
y0 = [-70.0;-70.0;0;0];
Tfin = 200;
h = 0.01;

t1 = 0:h:1;
n1 = length(t1)-1;
t2 = 1+h:h:Tfin-40;
n2 = length(t2)-1;
t3 = Tfin-40+h:h:Tfin;
n3 = length(t3)-1;
t = 0:h:Tfin;
newN = [n1 n2 n3];
% IAF params

Iext1_vec = [0.0 3.1 0.0]; %0.6;% 1.75;%1.13; %0.6; %1.1;%0.38; 
Iext2 = 0.0;

gc = 0.1;
gLexp = 0.0; %0.014; % for expIAF 0.014;
deltaT = 6.0; %[ 5.0 5.5 6.0];
V1 =[];
V2 = [];
G1_all =[];
G2_all =[];

% gL_vec = [0.01 0.05 0.1 0.5 0.8];
% gLexp = 0.5;
for j = 1:length(Iext1_vec)
     y = zeros(2,1);
     G1 = [];
     G2 = [];
     Iext1 = Iext1_vec(j);
     n = newN(j);
%     gLexp = gL_vec(j)
%     deltaT = deltaT_vec(j); %mV
    % parameters that get passed in:
    params = [C, gL, Eleak, vE, Iext1, Iext2, vT, gLexp, deltaT, gc];
    % this is to evaluate the derivative: dv/dt histV holds the sum over last
    % 50ms
    F = @(x,y,G) findDeriv_expIAF(x,y,G,params);
    if j == 1
        y(:,1) = y0(1:2);
        G1(1) = y0(3);
        G2(1) = y0(4);
    else
        y(:,1) = [V1(end);V2(end)];
        G1(1) = G1_all(end);
        G2(1) = G2_all(end);
    end
    for k = 2:n+1
        % conductance decays:
        G1(k) = G1(k-1)*exp(-h/sigmar);
        G2(k) = G2(k-1)*exp(-h/sigmar);
        G = [G1(k-1) G2(k-1)];
        
        %RK 2 step:
        K1 = h*F(t(k-1),y(:,k-1),G);  %K1 = h*f(tn-1,yn-1)
        K2 = h*F(t(k-1)+0.5*h, y(:,k-1) + (1/2)*K1, G);  %K2 = h*f(tn-1/2, yn + (1/2)*K1)
        K3 = h*F(t(k-1)+0.5*h, y(:,k-1) + (1/2)*K2, G); %K3 = h*f(tn+1/2, yn + (1/2)K2)
        K4 = h*F(t(k-1)+h, y(:,k-1) + K3,G);
        y(:,k) = y(:,k-1) + (1/6)*K1 + (1/3)*K2 + (1/3)*K3 + (1/6)*K4;
        
        % V1 spiked
        if  (y(1,k)>vT) %(y(1,k) > 80) 
            y(1,k) = vR;
            G2(k) = G2(k-1) + s;
        end
        % V1 spiked
        if  (y(2,k)>vT) %(y(2,k) > 80) 
            y(2,k) = vR;
            G1(k) = G1(k-1) + s;
        end    
    end
    V1 = [V1 y(1,:)];
    V2 = [V2 y(2,:)];
    
    G1_all = [G1_all G1];
    G2_all = [G2_all G2];
end
      
%%
figure
hold all
subplot(2,1,1)
hold all
plot(t,V1,'r','LineWidth',2.0)
axis([0,Tfin,-70,-45])
subplot(2,1,2)
hold all
plot(t,G2_all,'k','LineWidth',2.0)
axis([0,Tfin,0,0.1])

figure(2)
hold all
subplot(2,1,1)
hold all
plot(t,V2,'b','LineWidth',2.0)
axis([0,Tfin,-70,-45])
subplot(2,1,2)
hold all
plot(t,G1_all,'k','LineWidth',2.0)
axis([0,Tfin,0,0.1])


%% video -- just one neuron
figure
M = [];
dx = 100;
for i = 1:floor(length(V1)/dx)
   plot(t,-50*ones(size(t)),'--','color',[0.5 0.5 0.5],'LineWidth',2.0)
   hold on
     plot(t,vR*ones(size(t)),'--','color',[0.5 0.5 0.5],'LineWidth',2.0)
   hold on
   plot(t(1+(i-1)*dx:i*dx),V1(1+(i-1)*dx:i*dx),'k','LineWidth',3.0)
   hold on
   plot(t(i*dx),V1(i*dx),'o','MarkerSize',12,'MarkerFaceColor',[0.8 0.2 0.3],'Color',[0.1 0.1 0.1],'LineWidth',2.0)
   hold on
    title('Voltage')
    set(gca,'FontSize',25.0,'XAxisLocation','origin','YAxisLocation','origin','xTick',[-500 500],'yTick',[-500 500])
    axis([0,Tfin,-70,-45])
    M = [M,getframe(gcf)]; 
end

myVideo = VideoWriter(['oneNeuron_Exvoltage.avi']);
myVideo.FrameRate = 10;
open(myVideo);
writeVideo(myVideo, M);
close(myVideo);

%% video
figure
M = [];
dx = 100;
for i = 1:floor(length(V1)/dx)
   subplot(2,1,1)
   plot(t,-50*ones(size(t)),'k--','LineWidth',2.0)
   hold on
   plot(t(1+(i-1)*dx:i*dx),V1(1+(i-1)*dx:i*dx),'k','LineWidth',3.0)
   hold on
   plot(t(i*dx),V1(i*dx),'o','MarkerSize',12,'MarkerFaceColor',[0.8 0.2 0.3],'Color',[0.1 0.1 0.1],'LineWidth',2.0)
   hold on
    title('Voltage Neuron 1')
    set(gca,'FontSize',25.0,'XAxisLocation','origin','YAxisLocation','origin','xTick',[-500 500],'yTick',[-500 500])
    axis([0,Tfin,-70,-45])
    subplot(2,1,2)
        plot(t(1+(i-1)*dx:i*dx),V2(1+(i-1)*dx:i*dx),'k','LineWidth',3.0)
    hold on
    plot(t(i*dx),V2(i*dx),'o','MarkerSize',12,'MarkerFaceColor',[0.6 0.1 0.6],'Color',[0.1 0.1 0.1],'LineWidth',2.0)
    hold on
    title('Voltage Neuron 2')
    set(gca,'FontSize',25.0,'XAxisLocation','origin','YAxisLocation','origin','xTick',[-500 500],'yTick',[-500 500])
    axis([0,Tfin,-75,-55])  
    M = [M,getframe(gcf)]; 
end

myVideo = VideoWriter(['GJNeurons_voltage.avi']);
myVideo.FrameRate = 10;
open(myVideo);
writeVideo(myVideo, M);
close(myVideo);

   
