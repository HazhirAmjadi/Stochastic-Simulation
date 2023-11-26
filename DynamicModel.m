function [Out,Results,Max] = DynamicModel(StructureName,xis,mubar,omegas,inputname,plots)

global H

Stparameters = StructureModel(StructureName,xis,mubar,omegas);

switch inputname
    case 'ElCentro.xls'
        Xg=readtable(inputname);
        Xg = Xg.Acceleration;
    case 'WGN'
        global WGNX
        Xg = WGNX;
end

st0=0.02;
st=0.01;

T=(length(Xg)-1)*st0;
t=0:st:T;

switch StructureName
    case '5Story'
        I = 0.6;
    case '10Story'
        I = 1.2;
    case '20Story'
        I = 2;
end
xg(1:((T/st)+1))=0;
for i=1:length(xg)
    xg(i)=XG(Xg,st0,(i-1)*st)*9.807*I;
end
%% Results for system without TMD
Results.Ya = lsim(Stparameters.sysca,xg,t,Stparameters.q0);
Results.Yd = lsim(Stparameters.syscd,xg,t,Stparameters.q0);
Ydrift = zeros(size(Results.Yd,1),size(Results.Yd,2));
for i = 1 : size(Results.Yd,2)
    if i == 1 
        Ydrift(:,i) = Results.Yd(:,i);
    else
        Ydrift(:,i) = Results.Yd(:,i) - Results.Yd(:,i-1);
    end
end
Results.Ydrift = Ydrift./H;  

Results.Fd = Results.Yd.*Stparameters.k0./1000;
Shear= zeros(size(Results.Yd,1),size(Results.Yd,2));
for i = size(Results.Yd,2) : -1 : 1
    if i == size(Results.Yd,2) 
        Shear(:,i) = Results.Fd(:,i);
    else
        Shear(:,i) = Shear(:,i+1) + Results.Fd(:,i);
    end
end
Results.Shear = Shear;
%% Results for system with TMD

Results.Yta = lsim(Stparameters.syscta,xg,t,Stparameters.q0t);
Results.Ytd = lsim(Stparameters.sysctd,xg,t,Stparameters.q0t);
Ytdrift = zeros(size(Results.Yd,1),size(Results.Yd,2));
for i = 1 : size(Results.Yd,2)
    if i == 1 
        Ytdrift(:,i) = Results.Ytd(:,i);
    else
        Ytdrift(:,i) = Results.Ytd(:,i) - Results.Ytd(:,i-1);
    end
end
Results.Ytdrift = Ytdrift./H; 

Results.Ftd = Results.Ytd(:,1:size(Results.Ytd,2)-1).*Stparameters.k0./1000;
Sheart= zeros(size(Results.Yd,1),size(Results.Yd,2));
for i = size(Results.Yd,2) : -1 : 1
    if i == size(Results.Yd,2) 
        Sheart(:,i) = Results.Ftd(:,i);
    else
        Sheart(:,i) = Sheart(:,i+1) + Results.Ftd(:,i);
    end
end
Results.Sheart = Sheart;
%% Plot Results
PlotRe(Results,t,xg,StructureName,10,plots)
%% Get the maxResults.
% Max.woTMD(:,1): Maximum of the acceleration for structure without-tmd
% Max.woTMD(:,2): Maximum of the displacement for structure without-tmd
% Max.woTMD(:,3): Maximum of the drift for structure without-tmd
% Max.woTMD(:,4): Maximum of the Force for structure without-tmd
% Max.woTMD(:,5): Maximum of the Shear for structure without-tmd

% Max.TMD(:,1): Maximum of the acceleration for structure with-tmd
% Max.TMD(:,2): Maximum of the displacement for structure with-tmd
% Max.TMD(:,3): Maximum of the drift for structure with-tmd
% Max.TMD(:,4): Maximum of the Force for structure with-tmd
% Max.TMD(:,5): Maximum of the Shear for structure with-tmd
Max.woTMD = [max(abs(Results.Ya)./9.806)' max(abs(Results.Yd))' max(abs(Results.Ydrift))' max(abs(Results.Fd))' max(abs(Results.Shear))'];
Max.wTMD = [max(abs(Results.Yta)./9.806)' max(abs(Results.Ytd))' [max(abs(Results.Ytdrift)) 0]' [max(abs(Results.Ftd)) 0]' [max(abs(Results.Sheart)) 0]'];
% out = [TMD acc, TMD dis, TMD dif, ATMD/AwoTMD, disTMD/diswoTMD, difTMD/difwoTMD]
% Out = [max(Max.wTMD(:,1));max(Max.wTMD(:,2));max(Max.wTMD(:,3));max(Max.wTMD(:,1))/max(Max.woTMD(:,1));max(Max.wTMD(:,2))/max(Max.woTMD(:,2));max(Max.wTMD(:,3))/max(Max.woTMD(:,3))];
Out = max(Max.wTMD(:,1));
end

function g=XG(Xg,st,t)

%Xg= Acceleration vector
%st= sample time of accelerogram
%t=interest time
R=rem(t,st);
if R==0
    w=round(t/st)+1;
   g=Xg(w);
else
    mmm=round(((t-R)/st)+1);
    nnn=round(((t-R)/st)+2);
    a1=Xg(mmm);
    a2=Xg(nnn);
    g=R*((a2-a1)/(st))+a1;
end
end