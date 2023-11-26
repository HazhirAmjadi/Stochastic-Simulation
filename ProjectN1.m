clc
clear
close all
%%
rng(2)
global H 
StName = '5Story';                       % Structure name where show the number of storey that we consider for model
InputE = 'ElCentro.xls';                            % 'ElCentro.xls'

switch InputE
    case 'WGN'
    global WGNX
    WGNX = wgn(8001,1,1,1,2);
    WGNX = (WGNX./max(WGNX))*0.6;

end
H = 3.2;                                  % Height of each floor for calculate the failure denoted by drift
m_xis = 0.05;                             % The mean of TMD damping ratio
cv_xis = 0.01;                            % The sigma of TMD damping ratio
mu_xis = log((m_xis^2)/sqrt(cv_xis+m_xis^2));
sigma_xis = sqrt(log(cv_xis/(m_xis^2)+1));

m_mubar = 0.1;                           % The mean of TMD mass ratio mt/Ms, mt: tmd mass, Ms: structure mass
cv_mubar = 0.02;                         % The sigma of TMD damping ratio
mu_mubar = log((m_mubar^2)/sqrt(cv_mubar+m_mubar^2));
sigma_mubar = sqrt(log(cv_mubar/(m_mubar^2)+1));

m_ft = 1;                                % The mean of First mode frequency ratio
cv_ft = 0.05;                            % The sigma of First mode frequency ratio
mu_ft = log((m_ft^2)/sqrt(cv_ft+m_ft^2));
sigma_ft = sqrt(log(cv_ft/(m_ft^2)+1));

SimNumber = 1000;                        % Simulation number

%% Direct MCS
Lim = 0.63;                                % Acceleration Limit

Out = zeros(1,SimNumber);
If = zeros(1,SimNumber);
for i = 1 : SimNumber
    % Generate samples
    xis = lognrnd(mu_xis,sigma_xis);
    mubar = lognrnd(mu_mubar,sigma_mubar);
    ft = lognrnd(mu_ft,sigma_ft);
    Out(:,i) = DynamicModel(StName,xis,mubar,ft,InputE,'off');
    If(i) = Indicator(Out(:,i),Lim);
end
%% MCS Plot
% Estimate PF() by MCS disp(mean(Out'))
M1 = mean(Out);
CV1 = 1/sqrt(SimNumber)*std(Out)./mean(Out);
disp('Estimate mean by MCS for Acceleration:'); disp(M1)
disp('C.O.V for Acceleration:'); disp(CV1);
M1I = mean(If);
CV1I = 1/sqrt(SimNumber)*std(If)./mean(If);
disp('Estimate PF() by MCS for Acceleration:'); disp(M1I)
disp('C.O.V for Acceleration:'); disp(CV1I);

Iter = 1:SimNumber;
IMCS1 = cumsum(Out)./(Iter);
ICOVMCS1 = 1./sqrt(Iter).*sqrt(cumsum(Out.^2)./(Iter)-IMCS1.^2)./IMCS1;

IMCS1I = cumsum(If)./(Iter);
ICOVMCS1I = 1./sqrt(Iter).*sqrt(cumsum(If.^2)./(Iter)-IMCS1I.^2)./IMCS1I;

figure
subplot(2,1,1)
% plot(Iter,IMCS1)
% hold on
plot(Iter,IMCS1I)
xlabel('Iteration')
ylabel('Integral Estimation H')
axis tight
grid on
% legend('For Acceleration','For Failure')

subplot(2,1,2)
% % plot(Iter,ICOVMCS1)
% % hold on
plot(Iter,ICOVMCS1I)
axis tight
xlabel('Iteration')
ylabel('C.O.V Estimation')
grid on
% legend('For Acceleration','For Failure')
sgtitle('Based on MCS')
%% Convert now to standard gaussian space
DataNorm = randn(3,SimNumber);
Out2 = zeros(1,SimNumber);
If2 = zeros(1,SimNumber);

for i = 1 : SimNumber
    % Generate samples
    xis = exp(mu_xis+sigma_xis*DataNorm(1,i));
    mubar = exp(mu_mubar+sigma_mubar*DataNorm(2,i));
    ft = exp(mu_ft+sigma_ft*DataNorm(3,i));
    Out2(:,i) = DynamicModel(StName,xis,mubar,ft,InputE,'off');
    If2(i) = Indicator(Out2(:,i),Lim);
end
%% MCS Plot 2
% Estimate PF() by MCS disp(mean(Out'))
M2 = mean(Out2);
CV2 = 1/sqrt(SimNumber)*std(Out2)./mean(Out2);
disp('Estimate mean by MCS for Acceleration:'); disp(M2)
disp('C.O.V for Acceleration:'); disp(CV2);
M2I = mean(If2);
CV2I = 1/sqrt(SimNumber)*std(If2)./mean(If2);
disp('Estimate PF() by MCS for Acceleration:'); disp(M2I)
disp('C.O.V for PF() Acceleration:'); disp(CV2I);

Iter = 1:SimNumber;
IMCS2 = cumsum(Out2)./(Iter);
ICOVMCS2 = 1./sqrt(Iter).*sqrt(cumsum(Out2.^2)./(Iter)-IMCS2.^2)./IMCS2;

IMCS2I = cumsum(If2)./(Iter);
ICOVMCS2I = 1./sqrt(Iter).*sqrt(cumsum(If2.^2)./(Iter)-IMCS2I.^2)./IMCS2I;

figure
subplot(2,1,1)
% plot(Iter,IMCS2)
% hold on
plot(Iter,IMCS2I)
xlabel('Iteration')
ylabel('Integral Estimation H')
axis tight
grid on
% legend('For Acceleration','For Failure')

subplot(2,1,2)
% plot(Iter,ICOVMCS2)
% hold on
plot(Iter,ICOVMCS2I)
axis tight
xlabel('Iteration')
ylabel('C.O.V Estimation')
grid on
% legend('For Acceleration','For Failure')
sgtitle('Based on MCS with Standard Gaussian Space')

%% optimisation to get the maximum of the integrand
options = optimoptions(@fminunc,'Display','iter','PlotFcn','optimplotx','MaxFunctionEvaluations',1000);
input = fminunc(@(input) Optim(input,mu_xis,sigma_xis,mu_mubar,sigma_mubar,mu_ft,sigma_ft,StName,InputE),[0 0 0]',options);
xisOPT = exp(mu_xis+sigma_xis*input(1));
mubarOPT = exp(mu_mubar+sigma_mubar*input(2));
ftOPT = exp(mu_ft+sigma_ft*input(3));
disp('------------------------------------------')
disp('Optimom value for [Xis;Mu;Ft]:')
disp([xisOPT;mubarOPT;ftOPT])
%%
NoOpt = DynamicModel(StName,0.04,0.05,1,InputE,'anim');    % 'plot', 'anim' or off
Opt = DynamicModel(StName,xisOPT,mubarOPT,ftOPT,InputE,'anim');
%% IS method
mean_I=input;                    %From part C we have Importance Sampling mean
sigma_I=eye(3)*1;                  %From the solvings in the report

Out3 = zeros(1,SimNumber);
OutI3 = zeros(1,SimNumber);
If3 = zeros(1,SimNumber);
If3N = zeros(1,SimNumber);
for i = 1:SimNumber

    IS = mvnrnd(mean_I,sigma_I)';

    Q_I = mvnpdf(IS,[0;0;0],sigma_I)/mvnpdf(IS,mean_I,sigma_I);

    xisN = exp(mu_xis+sigma_xis*IS(1));
    mubarN = exp(mu_mubar+sigma_mubar*IS(2));
    ftN = exp(mu_ft+sigma_ft*IS(3));

    Out3(i) = DynamicModel(StName,xisN,mubarN,ftN,InputE,'off');
    If3(i) = Indicator(Out3(i),Lim);
    OutI3(i) = Q_I*Out3(i);
    If3N(i) = If3(i)*Q_I;
end

%% IS Plot
% Estimate PF() by MCS disp(mean(Out'))
M3 = mean(OutI3);
CV3 = 1/sqrt(SimNumber)*std(OutI3)./mean(OutI3);
disp('Estimate mean by IS for Acceleration:'); disp(M3)
disp('C.O.V for Acceleration by IS:'); disp(CV2);
M3I = mean(If3N);
CV3I = 1/sqrt(SimNumber)*std(If3N)./mean(If3N);
disp('Estimate PF() by IS for Acceleration:'); disp(M3I)
disp('C.O.V for Acceleration by IS:'); disp(CV3I);

IMCS3 = cumsum(OutI3)./(Iter);
ICOVMCS3 = 1./sqrt(Iter).*sqrt(cumsum(OutI3.^2)./(Iter)-IMCS3.^2)./IMCS3;

IMCS3I = cumsum(If3N)./(Iter);
ICOVMCS3I = 1./sqrt(Iter).*sqrt(cumsum(If3N.^2)./(Iter)-IMCS3I.^2)./IMCS3I;

figure
subplot(2,1,1)
plot(Iter,IMCS3)
hold on
plot(Iter,IMCS3I)
xlabel('Iteration')
ylabel('Integral Estimation H')
axis tight
grid on
legend('For Acceleration','For Failure')

subplot(2,1,2)
plot(Iter,ICOVMCS3)
hold on
plot(Iter,ICOVMCS3I)
axis tight
xlabel('Iteration')
ylabel('C.O.V Estimation')
grid on
legend('For Acceleration','For Failure')
sgtitle('Based on MCS')
%
figure
subplot(2,1,1)
plot(1:SimNumber,IMCS2I)
hold on
plot(1:SimNumber,IMCS3I)
ylabel('Integral Estimation H')
grid on
title('Comparison of MCS and IS')
legend('MCS','IS')
subplot(2,1,2)
plot(1:SimNumber,ICOVMCS2I)
hold on
plot(1:SimNumber,ICOVMCS3)
xlabel('Iteration (1:N)')
ylabel('C.O.V Estimation')
grid on
%% --------------------- [Indicator Function] ----------------------------

function If = Indicator(Out,Lim)

if Out(1) > Lim(1)
    If = 1;
else
    If = 0;
end
end


%% -------------------- [Optimization Function] ----------------------------

function Out = Optim(input,mu_xis,sigma_xis,mu_mubar,sigma_mubar,mu_ft,sigma_ft,StName,InputE)

xis = exp(mu_xis+sigma_xis*input(1));
mubar = exp(mu_mubar+sigma_mubar*input(2));
ft = exp(mu_ft+sigma_ft*input(3));
if mubar > 0.2
    Indx = 10;
else
    Indx = 1;
end
Outn = DynamicModel(StName,xis,mubar,ft,InputE,'off');
Out = Outn(1)*Indx;
end