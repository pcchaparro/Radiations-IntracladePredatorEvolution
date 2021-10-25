function dTr = EvoDyn(t,Tr,n,m,Food,Biom,Ea,landa,AMAX,THETAF,THETAP,tau,mu,sigma)

Trait=Tr;

Traitchange = zeros(m,1);
fGradF = zeros(m,1);
fGradP = zeros(m,1);
FitnessGrad = zeros(m,1);

for j=1:m
    %calculate fitness gradient (formula from symbolic solver
    %FitnessGrad.m)
    
    %From basal resources
    for i=1:n %number of resources
        fGradF(j,1) = fGradF(j,1) -(AMAX*Ea*Food(1,i)*exp(-(Trait(j,1) - THETAF(1,i))^2/(2*tau^2))*(2*Trait(j,1) - 2*THETAF(1,i)))/(2*tau^2);
    end
    %From feeding on other morphs (k)
    for k=1:m
        if k ~= j
            Biomdisp = landa*Biom(k,1);
            fGradP(j,1) = fGradP(j,1) -(AMAX*Ea*Biomdisp*exp(-(Trait(j,1) - THETAP)^2/(2*tau^2))*(2*Trait(j,1) - 2*THETAP))/(2*tau^2);
        end
    end
    
    FitnessGrad(j,1) = fGradF(j,1) + fGradP(j,1);
    Traitchange(j,1) = .5*mu*sigma*Biom(j,1)*FitnessGrad(j,1);
end

dTr = Traitchange;
end