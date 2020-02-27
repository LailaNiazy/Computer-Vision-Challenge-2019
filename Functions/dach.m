% berechnet die schiefsymmetrische Dachmatrix zum Vektor e
function e=dach(e)
e=[[0, -e(3), e(2)];[e(3),0,-e(1)];[-e(2),e(1),0]];
end