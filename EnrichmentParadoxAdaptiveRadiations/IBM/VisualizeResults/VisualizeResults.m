% This rutine load the data in the file 'Predator.txt' to produce the
% figure S5B. It can be used to read the outcome produced by the rutine
% IBMdiversification.

clear

% --- Figure S5B ---

data=importdata('Predator.txt');

n=10;
DTHETA=3;
TAU=1;
THETAP=0;
THETAF=(DTHETA*TAU)*(1:n);

histedges=-5.05:.1:max(THETAF)+5;
shist=size(histedges);
nbins=shist(1,2)-1;
sdata = size(data);
ecotrait=data(:,3:nbins+2);

figure
suptitle('Figure S5B')
title('Predator cannot evolve')
imagesc(1:sdata(1,1),histedges,flipud(log(ecotrait)'))
colormap(flipud(gray))
caxis([0 4])
yticks([THETAP THETAF])
ticklab = {'ThetaC'};
cont=0;
for j=n:-1:1
    cont = cont+1;
    ticklab{cont} = sprintf('Theta%d',j);
end
ticklab{cont+1} = 'ThetaC';
yticklabels(ticklab)
xlabel('Time')
ylabel('Niche Trait')