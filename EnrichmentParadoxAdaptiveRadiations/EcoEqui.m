function F = EcoEqui(t,x,AR,n,m,rho,FTmax,Ea,delta,landa)

    %Ecological dynamics
    % A is the matrix of attack rates of all ecomorphs on preexisting
    % resources and other ecomorphs. Rows correspond to each consumer ecomorph and
    % columns to the attacked resource. First column is the attack rate on
    % other ecomorphs
    % x is the vector of variables to solve
    
F = zeros(n+m,1);

Food=x(1:n,1);
Biom=x(n+1:n+m,1);

ingestF = zeros(m,1);
ingestP = zeros(m,1);
totingest= zeros(m,1);
biomcons= zeros(m,1);

Biomchange = zeros(m,1);
Foodchange = zeros(n,1);

for j=1:m %number of morphs
    %feeding on basal resources
    for i=1:n %number of resources
        ingestF(j,1) = ingestF(j,1) + AR(j,i+1)* Food(i,1); 
    end
    %feeding on other morphs (k)
    for k=1:m
        if k ~= j
            ingestP(j,1) = ingestP(j,1) + AR(j,1)*landa*Biom(k,1);
        end
    end
    totingest(j,1) = ingestF(j,1) + ingestP(j,1);
    %Biomass loss due to consumption by other morphs (k)
    for k=1:m
        if k ~= j
            biomcons(j,1)= biomcons(j,1) + AR(k,1)*landa*Biom(k,1);  
        end
    end
    Biomchange(j,1) = (Ea*totingest(j,1)-delta-biomcons(j,1))*Biom(j,1);
end

Fmax=FTmax/n;
grazing= zeros(n,1);
for i=1:n %number of food resources
    for j=1:m %number of morphs
        grazing(i,1) = grazing(i,1) + AR(j,i+1)*Food(i,1)*Biom(j,1);
    end
    Foodchange(i,1) = rho*(Fmax-Food(i,1)) - grazing(i,1);
end

F(1:n,1)=Foodchange;
F(n+1:n+m,1)=Biomchange;

end