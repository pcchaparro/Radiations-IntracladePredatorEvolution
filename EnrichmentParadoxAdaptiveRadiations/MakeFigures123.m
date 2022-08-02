% This script makes the figures 1, 2A and 3 of the paper "The enrichment
% paradox in adaptive radiations: emergence of predators hinders
% diversification in resource rich environments".
%
% Last modification: 25/10/21

clear

%%% --- Figure 1 ---

% Parameters
n=10;                     % Number of resources
rho=.01;                  % Resource growth rate
FTmax=10;                 % Total carrying capacity
AMAX=0.2;                 % Maximum attack rate
TAU=1;                    % Standard deviation of the Gaussian that determines how attack rate varies with trait
DTHETA=3;                 % Trait distance between optimals to feed on different resources
Ea=.6;                    % Efficiency to transform ingested mass into biomass
delta=.01;                % Mortality rate
landa=.1;                 % Proportion of biomass of the radation available for the intraclade consumer
mu = .001;                % mutation rate
sigma = .01;              % variance of mutational steps
fitland  = 1;             % if mutant fitness landscape will be display, it is necessary to specify the time at which the fitnes landscape is evaluated (before figures are generated)

% figure 1A (predator cannot evolve)

THETAP = -10; % Optimal trait value to feed upon other members of the radiation
RadiationFitnessLandscape

figure
suptitle('Figure 1')
tfig=0;
subplot(5,2,1)
for i=1:cont
    plot(tfig+(1:tspeciation{i}),TTrai{i},'k')
    hold on
    tfig = tfig + tspeciation{i};
end
xlim([0 1E4])
ylim([-2 32])
xlabel('Time')
ylabel('Niche trait')

% figure 1B (predator can evolve)

THETAP = 0; % Optimal trait value to feed upon other members of the radiation

% Run

RadiationFitnessLandscape

tfig=0;
subplot(5,2,2)
for i=1:cont
    plot(tfig+(1:tspeciation{i}),TTrai{i},'k')
    hold on
    tfig = tfig + tspeciation{i};
end
xlim([0 1E4])
ylim([-2 32])
xlabel('Time')
ylabel('Niche trait')

if fitland >0 && max(tfitland)<=timerad
    for i=1:stf(1,2)
        subplot(1+stf(1,2),2,[2*i+1 2*i+2])
        plot(mutrange,fitF(:,i))
        hold on
        plot([THETAP-5*TAU max(THETAF)+5*TAU],[0 0])
        xlim([-5 max(mutrange)])
        ylim([-.02 .14])
        xticks([THETAP THETAF])
        ticklab = {'ThetaC'};
        for j=1:n
            ticklab{j+1} = sprintf('Theta%d',j);
        end
        xticklabels(ticklab)
        yticklabels([0 .05 .1])
        title(sprintf('Time = %d',tfitland(1,i)))
        ylabel('Fitness')
    end
    xlabel('Niche trait')
end

%%% --- Figure 2A ---

tfig=0;
figure
suptitle('Figure 2A')
for i=1:cont
    h1=plot(tfig+(1:tspeciation{i}),sum(Predmort{i}),'k');
    hold on
    h2=plot(tfig+(1:tspeciation{i}),sum(Predmort{i})+delta,'r');
    tfig = tfig + tspeciation{i};
    legend([h1(1), h2(1)], 'Predator mortality','Total mortality')
end
xlim([0 1.5E4])
xlabel('Time')
ylabel('Mortality')

%%% --- Figure 3 ---

RadiationBiomDyn
     