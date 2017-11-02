clear variables

Nsim = 101;
phiMatr = zeros(Nsim,3);
extPhase = zeros(Nsim,1);
flexPhase = zeros(Nsim,1);

Ntot = 8;

Xold = rand(Ntot*2,1)';

da = linspace(0.005,0.9,Nsim);


% g(v) param
Vth = -50;
Vmax =  0;



Tfinale = 100000;

for k=1%:Nsim
    
    
    %networkParam
    Nactive = 8;
    Ntot = Nactive;
     
    gSynIn = 10;
    gSynEx = 10;
    
    EsynIn = -75;
    EsynEx = -10;
   
    tetaSyn = -30;%-31.6327;
    uSyn = 0.3;%0.2041
    
    syn = FTM_synapse_model(uSyn,tetaSyn);
    
    Di = zeros(Ntot,1);
    
    neuron = cell(Ntot,1);
    
    alpha =  da(k);
    
    % brainstream drive
    Di(1) = 0.1*alpha+0.0023;
    Di(2) = 0.1;
    Di(3) = 0.1*alpha+0.0023;
    Di(4) = 0.1;
    Di(5) = 0.104*alpha+0.001;
    Di(6) = 0.1;
    Di(7) = 0.104*alpha+0.001;
    Di(8) = 0.1;
    
    % Neuron with INaP parameters
    C = 10; %pF
    gNaP = 4.5; %nS
    ENa = 50; %mV
    gl = 4.5; %nS
    El = -62.5; %mV
    VhalfM = -40; %mV
    km = -6; %mV
    Vhalfh = -45; %mV
    kh = 4; %mV
    tau0 = 80; % ms
    tauMax = 160; % ms
    VhalfTau = -35; %mV
    kTau = 15; %mV
    
    % create neurons
    for i=1:Nactive
        neuron{i} = Danner_model(C,gNaP,ENa,gl,El,VhalfM,km,Vhalfh,kh,tau0,tauMax,VhalfTau,kTau,Di(i),gSynEx,EsynEx);
    end
    
    [g_in,g_ex,g_el] = setReducedConnectivityMatrix(Ntot,gSynIn,gSynEx,alpha);
    
    netw = CPG(Ntot,neuron,g_in,g_ex,g_el,EsynIn,EsynEx,syn,syn);
    
    
    
    % Simulate network
    
    
    
    opt.integrator = 'odeint';
    opt.Vth = -43;
    [T,X] = netw.sim([0 Tfinale],Xold,opt);
    phi = netw.getPhaseRepresentation([0 Tfinale],Xold,opt);
    TT{k} = T;
    XX{k} = X;
    Xold = X(end,:);
    try
    phiMatr(k,:) = phi(end,[2 4 6]);
    
    G = (X(:,1)-Vth)/(Vmax-Vth);
    G(G < 0) = 0;
    G(G > 1) = 1;
    
    iip = find(diff(G > 0.1) > 0,2,'last');
    iim = find(diff(G > 0.1) < 0,2,'last');
    
    if iim(end) > iip(end)
        flexPhase(k) = T(iim(1))-T(iip(1));
        extPhase(k) = T(iip(2))-T(iim(1));
    else
        flexPhase(k) = T(iim(2))-T(iip(1));
        extPhase(k) = T(iip(2))-T(iim(2));
    end
    catch
    end
        Tfinale = 5000;

    waitbar(k/Nsim);
    
end

close all

f = 1./(flexPhase+extPhase);

plot(phiMatr)


%%