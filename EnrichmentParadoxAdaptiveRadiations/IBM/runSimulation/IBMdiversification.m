% This rutine simulates a radiation using an individual-based model, with
% explicit genetic dynamics. Individuals are assigned two diploid
% genotypes: one set of loci determines the feeding niche trait and the
% other the mating trait. Mating is implemented following the approach of
% Dieckmann and Doebeli (1999). See supporting information of the paper
% "The enrichment paradox in adaptive radiations: emergence of predators
% hinders diversification in resource rich environments". 

clear

%Parameters----------------------------------------------------

n = 10;                    % Number of resources
rho=.01;                  % Resource growth rate
FTmax = 10;               % Total carrying capacity
AMAX=0.2;                 % Maximum attack rate
TAU=1;                    % Standard deviation of the Gaussian that determines how attack rate varies with trait
DTHETA=3;                 % Trait distance between optimals to feed on different resources
THETAP=0;                 % Optimal trait value to feed upon other members of the radiation
Ea=.6;                    % Efficiency to transform ingested mass into biomass
delta=.01;                % Mortality rate
landa=.1;                 % Proportion of biomass of the radation available for the intraclade consumer

nEloci = 10;              % number of loci of ecological trait
nAloci = 5;               % number of loci of assortative mating trait
mutrate = 0.005;          % mutation rate

volLake  = 1E3;           % Volumen of the habitat
tmax     = 1E6;           % Time points of the simulation
toutput  = 1E2;           % save population statistic every toutput time points
file_name= 'test.txt';    % file where statistics will be saved

%Initial conditions

niniind = 40;     %size of initial population
traitini= 4;      %initial trait of individuals

%Initialize food resources (assumed to be in its carrying capacity)

Fmaxi = FTmax/n*ones(1,n);
Fi = Fmaxi;
THETAF=(DTHETA*TAU)*(1:n);

%Create population matrix Pop: each row corresponds to each ind

Eloci = 1:2*nEloci;                   %columns 1:2*nEloci are ecological trait loci
Aloci = 2*nEloci+1:2*nEloci+2*nAloci; %columns 2*nEloci+1:2*nEloci+2*nAloci are assortative mating trait
sex   = 2*nEloci+2*nAloci+1;          %column  2*nEloci+2*nAloci+1 is sex
ecotr = 2*nEloci+2*nAloci+2;          %column  2*nEloci+2*nAloci+2 is ecological trait
mattr = 2*nEloci+2*nAloci+3;          %column  2*nEloci+2*nAloci+3 is assortative mating trait
Biores= 2*nEloci+2*nAloci+4;          %column  2*nEloci+2*nAloci+4 is reserve of biomass (when this reserve reaches 1 in a female, it produces 1 offspring)
aFi   = 2*nEloci+2*nAloci+5:2*nEloci+2*nAloci+5+n-1;          %other columns are the attack rate on resource ith
aP    = 2*nEloci+2*nAloci+5+n;        % column 2*nEloci+2*nAloci+5+n is the attack rate on other members of the radiation

nfemales = floor(niniind/2);
nmales   = niniind-nfemales;

Pop          = [ones(niniind,2*nEloci)*traitini/(2*nEloci) ones(niniind,nAloci) zeros(niniind,nAloci)]; % zeros(niniind,nAloci)]; %genotypes
Pop(:,sex)   = [ones(nfemales,1); -ones(nmales,1)]; %1 for females, and -1 for males
Pop(:,ecotr) = sum(Pop(:,Eloci),2);
Pop(:,mattr) = mean(Pop(:,Aloci),2);
Pop(:,Biores)= 0;
for j=1:n
    Pop(:,aFi(1,j)) = AMAX*exp(-((Pop(:,ecotr)-THETAF(1,j)).^2)./(2*TAU^2));
end
Pop(:,aP)    = AMAX*exp(-((Pop(:,ecotr)-THETAP).^2)./(2*TAU^2));

histedges=-5.05:.1:max(THETAF)+5;
shist=size(histedges);
nbins=shist(1,2)-1;

% create a file to save population

popfile = fopen(file_name,'w');

for i=1:tmax    
    
    %remove death individuals
    sPop = size(Pop);
    sp=sPop(1,1);
    Pop(rand(sp,1)<delta,:)=[];

    %population stats
    sPop = size(Pop);
    sp=sPop(1,1);
    nfemales = sum(Pop(:,sex)>0);
    nmales = sp-nfemales;
    
    %Feeding
    ingest = zeros(sp,n);
    for j=1:n
        ingest(:,j) = Pop(:,aFi(1,j)).*Fi(1,j); %Feeding on preexisting resources
    end
    ingestprey = Pop(:,aP).*landa*sp/volLake; %feeding on other individuals of the clade
    Pop(:,Biores)= Pop(:,Biores) + Ea.*(sum(ingest,2)+ingestprey); %energy available to use in biomass
    
    %update density of food resources
    for j=1:n
        Fi(1,j) = max(0,(Fi(1,j) + rho*(Fmaxi(1,j)-Fi(1,j)) - sum(ingest(:,j))/volLake));
    end
    
    %Reproduction
    Bioresfem = (Pop(:,sex)>0).*Pop(:,Biores);
    nrepfemale = sum(Bioresfem>=1);
    idrepfem=find(Bioresfem>=1); %id of reproducing females
    noffspring = floor(Bioresfem); %number of offspring per female
    Pop(:,Biores)= Bioresfem-noffspring; %remove the biomass used in reproduction
    idmales = find(Pop(:,sex)<0);
    Pop(idmales,Biores)= 0;      %reset males' reserve to avoid large numbers
    
    if nrepfemale>0
        %Create offspring matrix
        offspring = zeros(sum(noffspring),sPop(1,2));
        baby = 0;

        %find a dad for offspring

        ecotrmales = unique(Pop(idmales,ecotr));
        for j=1:nrepfemale
            mattrait = Pop(idrepfem(j,1),mattr);
            modmtr = mattrait*2-1; %normalized to compare with Dieckman&Doebeli,1999
            if modmtr>0
                sigmaA = 1/(20*modmtr^2); %from Dieckman&Doebeli,1999
                probdads = normpdf(ecotrmales,Pop(idrepfem(j,1),ecotr),sigmaA);
                ecotrsel = finddad(ecotrmales,probdads); %selected eco trait
                malesecotrsel = find((Pop(:,ecotr)>ecotrsel-1E-6 & Pop(:,ecotr)<ecotrsel+1E-6) & Pop(:,sex)<0);
                dad = randsample(malesecotrsel,1);
            elseif modmtr<0
                sigmaD = 1/(modmtr^2); %from Dieckman&Doebeli,1999
                probdads = 1-normpdf(ecotrmales,Pop(idrepfem(j,1),ecotr),sigmaD);
                ecotrsel = finddad(ecotrmales,probdads); %selected eco trait
                malesecotrsel = find((Pop(:,ecotr)>ecotrsel-1E-6 & Pop(:,ecotr)<ecotrsel+1E-6) & Pop(:,sex)<0);
                dad = randsample(malesecotrsel,1);
            else
                dad = randsample(idmales,1);
            end

            for k=1:noffspring(idrepfem(j,1),1)
                selalmom = unidrnd(2,1,nEloci+nAloci)+(0:2:(2*(nEloci+nAloci)-2)); %select genes from mom with full recombination
                selaldad = unidrnd(2,1,nEloci+nAloci)+(0:2:(2*(nEloci+nAloci)-2)); %select genes from dad with full recombination
                genfrommom = Pop(idrepfem(j,1),selalmom);
                genfromdad = Pop(dad,selaldad);
                genoff = [genfrommom;genfromdad];
                genoff = genoff(:)';
                baby = baby+1;
                offspring(baby,1:2*nEloci+2*nAloci) = genoff;
            end
        end

        %mutate the offspring randomly
        offspring(:,Aloci)=xor(offspring(:,Aloci),rand(sum(noffspring),2*nAloci)<mutrate);
        mutating=(rand(sum(noffspring),2*nEloci)<mutrate).*rand(sum(noffspring),2*nEloci);
        mutating(mutating<=.5&mutating>0)=-0.1;
        mutating(mutating>.5)=0.1;
        offspring(:,Eloci)=offspring(:,Eloci)+mutating;
   
        %complete the rest of the matrix
        offfemales = round(sum(noffspring)/2); %number of females in offspring
        offspring(:,sex) = [ones(offfemales,1); -ones(sum(noffspring)-offfemales,1)];
        offspring(:,ecotr) = sum(offspring(:,Eloci),2);
        offspring(:,mattr) = mean(offspring(:,Aloci),2);
        offspring(:,Biores)= 0;
        for j=1:n
            offspring(:,aFi(1,j)) = AMAX*exp(-((offspring(:,ecotr)-THETAF(1,j)).^2)./(2*TAU^2));
        end
        offspring(:,aP)    = AMAX*exp(-((offspring(:,ecotr)-THETAP).^2)./(2*TAU^2));

        %add the newborns to the population
        Pop = [Pop;offspring];
    end
    
    %Eliminate death individuals from predation
    predmort = sum(ingestprey);
    Pop(rand(sp,1)<(predmort/sp),:)=[];
    
    % write statistics every toutput time
    
    if rem(i,toutput)==0
        mattrbin = zeros(1,nbins);
        for j = 1:nbins
            idbin = find(Pop(:,ecotr)<histedges(1,j+1)&Pop(:,ecotr)>=histedges(1,j));
            if sum(idbin)>0
                mattrbin(1,j) = mean(Pop(idbin,mattr));
            end
        end
        dist=histcounts(Pop(:,ecotr),histedges); %count ind with trait ecotr in the bins of histedges
        sPop = size(Pop);
        sp=sPop(1,1);
        fprintf(popfile,'%.5g %.5g',i,sp);
        outputstr = ' %.5g';
        outputstr = repmat(outputstr, 1, 2*nbins);
        outputstr = [outputstr '\n'];
        fprintf(popfile,outputstr, dist, mattrbin);
    end
    
end