function fonctionInv=inverser_fonction(fonction)

fonctionInv=[];
for k=[1:length(fonction)]
    fonctionInv(end+1)=-fonction(k);
end
