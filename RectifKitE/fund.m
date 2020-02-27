function [F,el,er] = fund(pml,pmr)
%FUND Computes fundamental matrix and epipoles from camera matrices.
%
%[F,el,er] = fund(pml,pmr) calcola la matrice fondamentale
%F, l'epipolo sinistro el e destro er, partendo dalle due 
%matrici  di proiezione prospettica pml (MPP sinistra) e 
%pmr (MPP destra).

%calculates the fundamental matrix F, the left epipole el and the right
%epipole er, starting from the 2 projection matrices pml (left MPP) and pmr
%(right MPP)

%    Author: A. Fusiello 1999


%calcolo i centri ottici dalle due MPP
%Calculate the optical centers of both MPP
cl = -inv(pml(:,1:3))*pml(:,4);
cr = -inv(pmr(:,1:3))*pmr(:,4);

%calcolo gli epipoli come proiezione dei centri
%ottici
%calculate the epipoles with the projection of the optical centers
el = pml*[cr' 1]';
er = pmr*[cl' 1]';

%el = el./norm(el);
%er = er./norm(er);

%calcolo la matrice fondamentale
%calculate the fundamental matrix 
F=dach(er)*pmr(:,1:3)/(pml(:,1:3));

F = F./norm(F);

