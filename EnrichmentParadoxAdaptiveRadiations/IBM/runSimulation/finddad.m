function iddad=finddad(idmales,p)

%find n individuals randomly according to their probability p

pp=p/sum(p);
cp = [0, cumsum(pp)'];
ind = find(rand>cp, 1, 'last');
iddad = idmales(ind);

end