function Stparameters = StructureModel(StructureName,xis,mubar,omegas)

switch StructureName
    case '5Story'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        m0 = [12 12 12 11 10]*1e3;
        k0=[22e3 20e3 17.8e3 16e3 14.3e3]*1e3;
        xi=0.05;
        mv=1:length(m0);
        out1='acc';
        out2='dis';
        n=length(m0);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [ modeshape, Period, ~, ~, Freq] = ModalParameters( m0, k0 );
        [ mt,kt,ct] = TmdModal( m0,modeshape(:,1),xis,mubar,omegas*Freq(1));
        k=StiffnessMatrix(k0);
        c=DampingMatrix(diag(m0) ,k,xi);

        nd = length(m0)+1;

        Mn = zeros(nd);
        Mn(1:n,1:n) = diag(m0);
        Mn(nd,nd) = mt;

        Kn = zeros(nd);
        Kn(1:n,1:n) = k;
        Kn(nd,nd) = kt; Kn(n,n) = Kn(n,n)+kt; Kn(n,nd) = -kt; Kn(nd,n) = -kt;



        Cn = zeros(nd);
        Cn(1:n,1:n) = c;
        Cn(nd,nd) = ct;
        Cn(nd,nd) = ct; Cn(n,n) = Cn(n,n)+ct; Cn(n,nd) = -ct; Cn(nd,n) = -ct;
        mvt=1:length(Mn);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [Acta,Bcta,Cta,Dta]=StateSpaceAg(Mn,Kn,Cn,mvt,out1);
        syscta=ss(Acta,Bcta,Cta,Dta);
        q0t(1:2*nd)=0;
        q0t=q0t';

        [Actd,Bctd,Ctd,Dtd]=StateSpaceAg(Mn,Kn,Cn,mvt,out2);
        sysctd=ss(Actd,Bctd,Ctd,Dtd);

        [Aca,Bca,Ca,Da]=StateSpaceAg(diag(m0),k,c,mv,out1);
        sysca=ss(Aca,Bca,Ca,Da);

        [Acd,Bcd,Cd,Dd]=StateSpaceAg(diag(m0),k,c,mv,out2);
        syscd=ss(Acd,Bcd,Cd,Dd);

        q0(1:2*n)=0;
        q0=q0';

    case '10Story'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        m0 = [179 170 161 152 143 134 125 116 107 98]*1e3;
        k0=[62.47 59.26 56.14 53.02 49.91 46.79 43.67 40.55 37.43 34.31]*1e6;
        xi=0.05;
        mv=1:length(m0);
        out1='acc';
        out2='dis';
        n=length(m0);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [ modeshape, Period, ~, ~, Freq] = ModalParameters( m0, k0 );
        [ mt,kt,ct] = TmdModal( m0,modeshape(:,1),xis,mubar,omegas*Freq(1));
        k=StiffnessMatrix(k0);
        c=DampingMatrix(diag(m0) ,k,xi);

        nd = length(m0)+1;

        Mn = zeros(nd);
        Mn(1:n,1:n) = diag(m0);
        Mn(nd,nd) = mt;

        Kn = zeros(nd);
        Kn(1:n,1:n) = k;
        Kn(nd,nd) = kt; Kn(n,n) = Kn(n,n)+kt; Kn(n,nd) = -kt; Kn(nd,n) = -kt;



        Cn = zeros(nd);
        Cn(1:n,1:n) = c;
        Cn(nd,nd) = ct;
        Cn(nd,nd) = ct; Cn(n,n) = Cn(n,n)+ct; Cn(n,nd) = -ct; Cn(nd,n) = -ct;
        mvt=1:length(Mn);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [Acta,Bcta,Cta,Dta]=StateSpaceAg(Mn,Kn,Cn,mvt,out1);
        syscta=ss(Acta,Bcta,Cta,Dta);
        q0t(1:2*nd)=0;
        q0t=q0t';

        [Actd,Bctd,Ctd,Dtd]=StateSpaceAg(Mn,Kn,Cn,mvt,out2);
        sysctd=ss(Actd,Bctd,Ctd,Dtd);

        [Aca,Bca,Ca,Da]=StateSpaceAg(diag(m0),k,c,mv,out1);
        sysca=ss(Aca,Bca,Ca,Da);

        [Acd,Bcd,Cd,Dd]=StateSpaceAg(diag(m0),k,c,mv,out2);
        syscd=ss(Acd,Bcd,Cd,Dd);

        q0(1:2*n)=0;
        q0=q0';

    case '20Story'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % m0 = [283 276 276 276 276 276 276 276 276 276 276 276 276 276 276 276 276 276 276 276]*1e3;
        % k0=[136.92 136.92 135.94 135.83 134.67 133.38 133.07 132.72 132.31 131.83 125.29 118.13 116.06 108.71 103.37 100.24 91.88 80.68 64.86 45.94]*1e6;
        m0 = [140 140 140 140 140 140 140 140 140 140]*1e3;
        k0=[220.6 220.6 220.6 192.5 192.5 192.5 176.1 176.1 176.1 164.3]*1e6;
        xi=0.05;
        mv=1:length(m0);
        out1='acc';
        out2='dis';
        n=length(m0);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [ modeshape, Period, ~, ~, Freq] = ModalParameters( m0, k0 );
        [ mt,kt,ct] = TmdModal( m0,modeshape(:,1),xis,mubar,omegas*Freq(1));
        k=StiffnessMatrix(k0);
        c=DampingMatrix(diag(m0) ,k,xi);

        nd = length(m0)+1;

        Mn = zeros(nd);
        Mn(1:n,1:n) = diag(m0);
        Mn(nd,nd) = mt;

        Kn = zeros(nd);
        Kn(1:n,1:n) = k;
        Kn(nd,nd) = kt; Kn(n,n) = Kn(n,n)+kt; Kn(n,nd) = -kt; Kn(nd,n) = -kt;



        Cn = zeros(nd);
        Cn(1:n,1:n) = c;
        Cn(nd,nd) = ct;
        Cn(nd,nd) = ct; Cn(n,n) = Cn(n,n)+ct; Cn(n,nd) = -ct; Cn(nd,n) = -ct;
        mvt=1:length(Mn);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [Acta,Bcta,Cta,Dta]=StateSpaceAg(Mn,Kn,Cn,mvt,out1);
        syscta=ss(Acta,Bcta,Cta,Dta);
        q0t(1:2*nd)=0;
        q0t=q0t';

        [Actd,Bctd,Ctd,Dtd]=StateSpaceAg(Mn,Kn,Cn,mvt,out2);
        sysctd=ss(Actd,Bctd,Ctd,Dtd);

        [Aca,Bca,Ca,Da]=StateSpaceAg(diag(m0),k,c,mv,out1);
        sysca=ss(Aca,Bca,Ca,Da);

        [Acd,Bcd,Cd,Dd]=StateSpaceAg(diag(m0),k,c,mv,out2);
        syscd=ss(Acd,Bcd,Cd,Dd);

        q0(1:2*n)=0;
        q0=q0';
end

Stparameters.q0t = q0t;
Stparameters.q0 = q0;
Stparameters.syscta = syscta;
Stparameters.sysctd = sysctd;
Stparameters.sysca = sysca;
Stparameters.syscd = syscd;
Stparameters.modeshape = modeshape;
Stparameters.Period = Period;
Stparameters.Freq = Freq;
Stparameters.k0 = k0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c=DampingMatrix(m ,k,xi)

[a, b]=eig(m\k);
bb=sqrt(b);
[ phi, omega ] = sortmode( a, bb );
M=phi'*m*phi;
cc=2*xi*omega*M;
c=phi'\cc/phi;

end

function [ modeshape, Period, MP, MPsum, omega,Frequency] = ModalParameters( m0, k0 )
%m : Mass vector
%k: stiffness vector

n=length(m0);
L=ones(n,1);
m=diag(m0);
k=zeros(n,n);
for i=1:n

    if i<n
        k(i,i)=k0(i)+k0(i+1);
        k(i,i+1)=-k0(i+1);
        k(i+1,i)=-k0(i+1);
    else
        k(i,i)=k0(i);
    end
end
[fi, om]=eig(m\k);
[ phi, omega ] = sortmode( fi, sqrt(om) );
Frequency=diag(omega)/(2*pi);
Period=1./Frequency;
for i=1:n
    phi(:,i)= phi(:,i)/phi(n,i);
end
modeshape=phi;
M=phi'*m*phi;
eta=phi'*m*L;
mp=(eta.^2)./diag(M);
MP=100*mp/sum(m0);
MPsum=MP;
for i=1:n
    if i>1
        MPsum(i)=sum(MP(1:i));
    end

end

end

function [ modeshape, frequency ] = sortmode( phi, omega )

n=length(diag(omega));
modeshape=zeros(n,n);
frequency=zeros(n,n);
omegav=diag(omega);
for i=1:n
    [frequency(i,i), z]=min(omegav);
    modeshape(:,i)=phi(:,z);
    sum=0;
    for j=1:length(omegav)
        if j==z
        else
            sum=sum+1;
            omegad(sum)=omegav(j);
            phid(:,sum)=phi(:,j);
        end
    end
    clear omegav phi
    omegav=omegad;
    phi=phid;
    if length(omegav)>1
        clear omegad phid
    end
end
end

function [A,B,C,D]=StateSpaceAg(m,k,c,mv,out)
%systemmatrices function produce state space matrices fo continues system
%m= Mass Matrix
%k= Stiffness Matrix
%c= Damping Matrix
%mv= location vector of output


mv=sort(mv);
mo=length(mv);


n=length(m);



A(1:n,1:n)=zeros(n,n);
A(1:n,n+1:2*n)=eye(n);
A(n+1:2*n,1:n)=-inv(m)*k;
A(n+1:2*n,n+1:2*n)=-inv(m)*c;

L=ones(n,1);


B(1:n,1)=0;
B(n+1:2*n,1)=-L;


In=eye(n);
Cs(1:mo,1:n)=0;
for i=1:mo
    Cs(i,:)=In(mv(i),:);
end

C(1:mo,1:2*n)=0;
if strcmp(out,'dis')==1
    C(1:mo,1:n)=Cs;
    D=zeros(mo,1);
end
if strcmp(out,'vel')==1
    C(1:mo,n+1:2*n)=Cs;
    D=zeros(mo,1);
end
if strcmp(out,'acc')==1
    C(1:mo,1:n)=-Cs*inv(m)*k;
    C(1:mo,n+1:2*n)=-Cs*inv(m)*c;
    D(1:mo,1)=-Cs*L;
end

end

function k=StiffnessMatrix(k0)
n=length(k0);
k=zeros(n,n);
for i=1:n
    if i<n
        k(i,i)=k0(i)+k0(i+1);
        k(i,i+1)=-k0(i+1);
        k(i+1,i)=-k0(i+1);
    else
        k(i,i)=k0(i);
    end
end

end
function [ mt,kt,ct] = TmdModal( ms,phi1,xis,mubar,omegas)

n=length(phi1);
L=ones(n,1);
phibar=phi1*(phi1'*diag(ms)*L)/(phi1'*diag(ms)*phi1);
Mbar1=phibar'*diag(ms)*phibar;
phibarn=phibar(n);
f=(1/(1+mubar*phibarn))*(1-xis*sqrt(mubar*phibarn/(1+mubar*phibarn)));
xit=phibarn*((xis/(1+mubar))+(sqrt(mubar/(1+mubar))));
mt=mubar*Mbar1;
kt=(f^2)*(omegas^2)*mt;
ct=2*xit*mt*sqrt(kt/mt);

end

