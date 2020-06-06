function [] = elementcount(Tnod, nelements)
%%% ESTO SIRVE PARA CONTAR SI HAY ALGÚN ELEMENTO REPETIDO

elements_disc = zeros(nelements, 3);
for i=1:length(Tnod)
    n1 = Tnod(i,1);
    n2 = Tnod(i,2);
    found = 0;
    for k = 1:length(elements_disc)
            while found == 0 
                if (n1 == elements_disc(i,1) && n2 == elements_disc(i,2)) || (n2 == elements_disc(i,1) && n1 == elements_disc(i,2))
                    elements_disc(i,3) = elements_disc(i,3)+1;
                    found = 1;
                    disp('Hay un elemento repetido');
                else
                    elements_disc(i,1) = Tnod(i,1);
                    elements_disc(i,2) = Tnod(i,2);
                    elements_disc(i,3) = 1;
                    found = 1;
                end
            end
    end
end

a = 1;
end