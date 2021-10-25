% This rutine computes the temporal dynamics of a radiation seeded with
% an ancestal population.

%--------------------------------------------------------------------
%ANCESTRAL POPULATION

%Initial conditions
initialFood = FTmax/n*ones(n,1);
initialBiom = .1;
initialTrai = 2;

THETAF=(DTHETA*TAU)*(1:n); %Optimal trait values for n resources

m=1;
tmax0 = 5000;

%eco-evo dynamics

tmax1 = 10000;
Food = zeros(n,tmax1);
Biom = zeros(m,tmax1);
Trai = zeros(m,tmax1);

Food(:,1) = initialFood;
Biom(:,1) = initialBiom;
Trai(:,1) = initialTrai;

AR   = zeros(m,n+1);
traitchange = 1;
time=1;

while (time<tmax1-1 && traitchange>1E-8)
    for j=1:m
    %Attack rates
        AR(j,1) = AMAX.*exp(-((Trai(j,time)-THETAP).^2)./(2*(TAU^2))); %attack rate on other ecomorphs
        AR(j,2:end) = AMAX.*exp(-((Trai(j,time)-THETAF).^2)./(2*(TAU^2))); %attack rate of morph j on preexisting resources
    end
    
    %find the ecological equilibrium
    tmaxe = 100;
    dif = 1;
    rep = 0;
    x0    = [Food(:,time)' Biom(:,time)'];
    while (dif>1E-6 && rep<20)
        fod   = @(t,x) EcoEqui(t,x,AR,n,m,rho,FTmax,Ea,delta,landa);
        [t1,x1] = ode45(fod,[0 tmaxe],x0);
        dif = abs(x1(end-1,1)-x1(end,1));
        rep = rep + 1;
        x0  = x1(end,:);
    end

    Food(:,time+1)=x1(end,1:n);
    Biom(:,time+1)=x1(end,n+1:n+m);
    
    fod   = @(t,Tr) EvoDyn(t,Tr,n,m,Food(:,time+1)',Biom(:,time+1),Ea,landa,AMAX,THETAF,THETAP,TAU,mu,sigma);
    [tevo,xevo] = ode23(fod,[0 1E6],Trai(:,time));
    Trai(:,time+1) = xevo(end,:);
    
    traitvar = zeros(m,1);
    for j=1:m
        traitvar(j,1) = abs(Trai(j,time+1)-Trai(j,time));
    end
    traitchange = max(traitvar);
    time = time+1;
end

TFood{1} = Food(:,1:time);
TBiom{1} = Biom(:,1:time);
TTrai{1} = Trai(:,1:time);
Predmort{1} = zeros(1,time);
tspeciation{1} = time;

%Calculate the 2nd derivative to determine if selection is disruptive

disruptive = zeros(1,m);
for j=1:m
    secderF = 0;
    %from preexisting resources
    for i=1:n
        secderF = secderF -(AMAX*Ea*TFood{1}(i,end)*exp(-(TTrai{1}(j,end) - THETAF(1,i))^2/(2*TAU^2))*(- TTrai{1}(j,end)^2 + 2*TTrai{1}(j,end)*THETAF(1,i) + TAU^2 - THETAF(1,i)^2))/TAU^4;
    end
    %From feeding on other morphs (k)
    secderP = 0;
    for k=1:m
        if k ~= j
            Biomdisp = landa*TBiom{1}(k,end);
            secderP = secderP -(AMAX*Ea*Biomdisp*exp(-(TTrai{1}(j,end) - THETAP)^2/(2*TAU^2))*(- TTrai{1}(j,end)^2 + 2*TTrai{1}(j,end)*THETAP + TAU^2 - THETAP^2))/TAU^4;
        end
    end
    secder = secderF + secderP;
    if secder>0
        disruptive(1,j)=1;
    else
        disruptive(1,j)=0;
    end
end

%--------------------------------------------------------------------
%SPECIATION
%If selection is disruptive split the population(s) with disruptive
%selection

timemaxrad=1E8;
timerad=0;
cont = 1;
numsp = 1;
numspnopred = 1;

while (sum(disruptive)>0 && m<n+15) && timerad<timemaxrad
    
    timerad = timerad + time;
    
    initialFood = TFood{cont}(:,end);
    initialBiom = [];
    initialTrai = [];
    
    %Diversification

    for j=1:m
        if disruptive(1,j)==1
            initialBiom = [initialBiom; TBiom{cont}(j,end)/2; TBiom{cont}(j,end)/2];
            initialTrai = [initialTrai; TTrai{cont}(j,end)*1-1E-3; TTrai{cont}(j,end)*1+1E-3];
        else
            initialBiom = [initialBiom; TBiom{cont}(j,end)];
            initialTrai = [initialTrai; TTrai{cont}(j,end)];
        end
    end
    
    m=m+sum(disruptive);
    
    % Eco-evo dynamics
    
    tmax1 = 10000;
    Food = zeros(n,tmax1);
    Biom = zeros(m,tmax1);
    Trai = zeros(m,tmax1);

    Food(:,1) = initialFood;
    Biom(:,1) = initialBiom;
    Trai(:,1) = initialTrai;

    AR   = zeros(m,n+1);
    traitchange = 1;
    time=1;
    
    while (time<tmax1-1 && traitchange>1E-8)
        for j=1:m
        %Attack rates
            AR(j,1) = AMAX.*exp(-((Trai(j,time)-THETAP).^2)./(2*(TAU^2))); %attack rate on other ecomorphs
            AR(j,2:end) = AMAX.*exp(-((Trai(j,time)-THETAF).^2)./(2*(TAU^2))); %attack rate of morph j on preexisting resources
        end

        %find the ecological equilibrium
        tmaxe = 100;
        dif = 1;
        rep = 0;
        x0    = [Food(:,time)' Biom(:,time)'];
        while (dif>1E-6 && rep<20)
            fod   = @(t,x) EcoEqui(t,x,AR,n,m,rho,FTmax,Ea,delta,landa);
            [t1,x1] = ode45(fod,[0 tmaxe],x0);
            dif = abs(x1(end-1,1)-x1(end,1));
            rep = rep + 1;
            x0  = x1(end,:);
        end

        Food(:,time+1)=x1(end,1:n);
        Biom(:,time+1)=x1(end,n+1:n+m);

        fod   = @(t,Tr) EvoDyn(t,Tr,n,m,Food(:,time+1)',Biom(:,time+1),Ea,landa,AMAX,THETAF,THETAP,TAU,mu,sigma);
        [tevo,xevo] = ode23(fod,[0 1E6],Trai(:,time));
        Trai(:,time+1) = xevo(end,:);

        traitvar = zeros(m,1);
        for j=1:m
            traitvar(j,1) = abs(Trai(j,time+1)-Trai(j,time));
        end
        traitchange = max(traitvar);
        time = time+1;
    end

    TFood{cont+1} = Food(:,1:time);
    TBiom{cont+1} = Biom(:,1:time);
    TTrai{cont+1} = Trai(:,1:time);
    tspeciation{cont+1} = time;
    
    %Calculate predation mortality
    
    for j=1:m
        ARsim = AMAX.*exp(-((TTrai{cont+1}-THETAP).^2)./(2*(TAU^2))); %attack rate on other ecomorphs
        Predmort{cont+1} = landa.*ARsim.*TBiom{cont+1};
    end

    %Calculate the 2nd derivative to determine if selection is disruptive

    disruptive = zeros(1,m);
    for j=1:m
        secderF = 0;
        %from preexisting resources
        for i=1:n
            secderF = secderF -(AMAX*Ea*TFood{cont+1}(i,end)*exp(-(TTrai{cont+1}(j,end) - THETAF(1,i))^2/(2*TAU^2))*(- TTrai{cont+1}(j,end)^2 + 2*TTrai{cont+1}(j,end)*THETAF(1,i) + TAU^2 - THETAF(1,i)^2))/TAU^4;
        end
        %From feeding on other morphs (k)
        secderP = 0;
        for k=1:m
            if k ~= j
                Biomdisp = landa*TBiom{cont+1}(k,end);
                secderP = secderP -(AMAX*Ea*Biomdisp*exp(-(TTrai{cont+1}(j,end) - THETAP)^2/(2*TAU^2))*(- TTrai{cont+1}(j,end)^2 + 2*TTrai{cont+1}(j,end)*THETAP + TAU^2 - THETAP^2))/TAU^4;
            end
        end
        secder = secderF + secderP;
        if secder>0
            disruptive(1,j)=1;
        else
            disruptive(1,j)=0;
        end
    end
    cont = cont+1;
end

timerad = timerad + time;

%------------------------------------------------------------
% MUTANT FITNESS LANDSCAPE

%times at which the fitness landscape should be evaluated:
tfitland = [tspeciation{1} 2600 6120 timerad];

%Find the points in the simulation
stf = size(tfitland);
ite = zeros(stf(1,2),1);
tite = zeros(stf(1,2),1);
mutrange = THETAP-5*TAU:.1:max(THETAF)+5*TAU; %trait values to be evaluated
smr = size(mutrange);
% nopred = 1:smr(1,2);
nopred = find(mutrange>THETAP+2*TAU);
spred = size(nopred);

if fitland >0 && max(tfitland)<=timerad
    for j=1:stf(1,2)
        cti = 0;
        foal = 0;
        for i=1:cont
            cti = cti + tspeciation{i};
            if tfitland(1,j)<=cti && foal == 0
                ite(j,1) = i;
                tite (j,1) = tfitland(1,j) - (cti - tspeciation{i});
                foal = 1;
            end
        end
    end
    
    fitF = zeros(smr(1,2),stf(1,2));
    for k=1:stf(1,2)
        sTB = size(TBiom{ite(k,1)}(:,tite(k,1)));
        numec = sTB(1,1);
        ingest = zeros(1,smr(1,2));
        %feeding on basal resources
        for i=1:n %number of resources
            ARFood = AMAX*exp(-(mutrange-THETAF(1,i)).^2/(2*(TAU^2)));
            ingest = ingest + ARFood*TFood{ite(k,1)}(i,tite(k,1));
        end
        mort = delta*ones(1,smr(1,2));
        
        %feeding on other ecomorphs
        if numec>1
            ARPrey = AMAX*exp(-(mutrange-THETAP).^2/(2*(TAU^2)));
            Biomav = sum(TBiom{ite(k,1)}(2:end,tite(k,1)));
            ingest = ingest + ARPrey*landa*Biomav;
        
        %loss due to predation mortality      
            for l=1:numec
                ARPred = AMAX*exp(-(TTrai{ite(k,1)}(l,tite(k,1))-THETAP).^2/(2*TAU^2));
                mort(nopred) = mort(nopred) + ARPred*landa*TBiom{ite(k,1)}(l,tite(k,1));
            end
        end
        fitF(:,k) = Ea*ingest-mort;
    end
end