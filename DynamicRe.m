function [Results,Max] = DynamicRe(StructureName,xis,mubar,omegas,inputname,plots)

Stparameters = StructureModel(StructureName,xis,mubar,omegas);

Xg=readtable(inputname);
Xg = Xg.Acceleration;
st0=0.02;
st=0.01;

T=(length(Xg)-1)*st0;
t=0:st:T;


xg(1:((T/st)+1))=0;
for i=1:length(xg)
    xg(i)=XG(Xg,st0,(i-1)*st)*9.807;
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
Results.Ydrift = Ydrift; 

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
Results.Ytdrift = Ytdrift; 

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
%% Get the maxResults
Max.woTMD = [max(abs(Results.Ya./9.806))' max(abs(Results.Yd))' max(abs(Results.Ydrift))' max(abs(Results.Fd))' max(abs(Results.Shear))'];
Max.wTMD = [max(abs(Results.Yta./9.806))' max(abs(Results.Ytd))' [max(abs(Results.Ytdrift)) 0]' [max(abs(Results.Ftd)) 0]' [max(abs(Results.Sheart)) 0]'];
%% Plot Results
PlotRe(Results,t,xg,StructureName,10,plots)
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