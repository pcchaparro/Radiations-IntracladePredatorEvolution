clear

%Parameters----------------------------------------------------

n =10;                    % Number of resources
rho=.01;                  % Resource growth rate
FTmax = 10;               % Total carrying capacity
AMAX=0.2;                 % Maximum attack rate
TAU=1;                    % Standard deviation of the Gaussian that determines how attack rate varies with trait
DTHETA=3;                 % Trait distance between optimals to feed on different resources
THETAP=[-10 0];           % Optimal trait value to feed upon other members of the radiation
Ea=.6;                    % Efficiency to transform ingested mass into biomass
delta=.01;                % Mortality rate
landa=.1;                 % Proportion of biomass of the radation available for the intraclade consumer
mu = .001;                % mutation rate
sigma = .01;              % variance of mutational steps

FTmaxV  = [5 8 10];
THETAPV = [-10 0];

sFT = size(FTmaxV);
sTP = size(THETAPV);
ite = 1;

for jF=1:sTP(1,2)
    THETAP=THETAPV(1,jF);
    
    for iF=1:sFT(1,2)
        FTmax = FTmaxV(1,iF);
        
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

        TFood{1} = Food(:,1:time);      % Resource density
        TBiom{1} = Biom(:,1:time);      % Biomass of ecomorphs
        TTrai{1} = Trai(:,1:time);      % Traits of ecomorphs
        TNoPr{1} = Biom(:,1:time);      % Biomass of non-predators
        TnNoP{1} = ones(1,time);        % Number of non-predator ecomorphs
        TPred{1} = zeros(1,time);       % Biomass of predators
        TnPre{1} = zeros(1,time);       % Number of predator ecomorphs
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

        timemaxrad=3E4;
        timerad=0;
        cont = 1;
        numsp = 1;
        numspnopred = 1;

        while (sum(disruptive)>0 && m<n+10) && timerad<timemaxrad

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
            TNoPr{cont+1} = [];
            TPred{cont+1} = [];
            TnNoP{cont+1} = [];
            tspeciation{cont+1} = time;

            pred = 0;
            for j=1:m
                if Trai(j,time) > THETAP+.5*DTHETA
                    TNoPr{cont+1} = [TNoPr{cont+1}; Biom(j,1:time)];
                    TnNoP{cont+1} = [TnNoP{cont+1}; ones(1,time)];
                else
                    TPred{cont+1} = [TPred{cont+1}; Biom(j,1:time)];
                    pred = 1;
                end
            end

            if pred >0
                TnPre{cont+1} = ones(1,time);
            else
                TnPre{cont+1} = zeros(1,time);
                TPred{cont+1} = zeros(1,time);
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

        tfig=0;
        tt=[];
        sp=[];
        sn=[];
        nn=[];
        np=[];
        for i=1:cont
            sTPred = size(TPred{i});
            sTNoPr = size(TNoPr{i});
            if sTPred(1,1) == 1
                sp = [sp TPred{i}];
            else
                sp = [sp sum(TPred{i})];
            end
            if sTNoPr(1,1) == 1
                sn = [sn TNoPr{i}];
                nn = [nn TnNoP{i}];
            else
                sn = [sn sum(TNoPr{i})];
                nn = [nn sum(TnNoP{i})];
            end
            np = [np TnPre{i}];
            tt = [tt tfig+(1:tspeciation{i})];
            tfig = tfig + tspeciation{i};
        end
        BiomPred{ite} = sp;
        BiomNoPred{ite} = sn;
        Timerad{ite} = tt;
        nPred{ite} = np;
        nNoPred{ite} = nn;
        ite = ite+1;
    end
end

figure
suptitle('Figure 3')

for j = 1:sTP(1,2)
    
    subplot(1+sTP(1,2),2,j)
    for i=1:sFT(1,2)
        hold on
        plot(Timerad{(j-1)*sFT(1,2)+i},nNoPred{(j-1)*sFT(1,2)+i},'linewidth',2) 
        plot(Timerad{(j-1)*sFT(1,2)+i},nPred{(j-1)*sFT(1,2)+i},'linewidth',2)
    end
    ylim([0 11])
    xlabel('Time')
    ylabel('Number of ecomorphs')
    
    subplot(1+sTP(1,2),2,sTP(1,2)+j)
    for i=1:sFT(1,2)
        plot(Timerad{(j-1)*sFT(1,2)+i},BiomNoPred{(j-1)*sFT(1,2)+i},'linewidth',2)
        hold on
        plot(Timerad{(j-1)*sFT(1,2)+i},BiomPred{(j-1)*sFT(1,2)+i},'linewidth',2)   
    end
    ylim([0 6])
    xlabel('Time')
    ylabel('Abundance')
end

%subplot 5-6 are only to enable the visualization of the legend
h2=subplot(1+sTP(1,2),2,[5 6]);
for i=1:sFT(1,2)
    hold on
    plot(h2,Timerad{(j-1)*sFT(1,2)+i},nNoPred{(j-1)*sFT(1,2)+i},'linewidth',2,'Visible','off') 
    plot(h2,Timerad{(j-1)*sFT(1,2)+i},nPred{(j-1)*sFT(1,2)+i},'linewidth',2,'Visible','off')
end
set(h2,'Visible','off');
leg1={'Low productivity - other ecomorphs','Low productivity - intraclade consumer','Intermediate productivity - other ecomorphs','Intermediate productivity - intraclade consumer','High productivity - other ecomorphs','High productivity - intraclade consumer'};
lgd1=legend (leg1,'location','southoutside');