# Computer-Vision-Challenge-2019


Die Computer Vision Challenge besteht in diesem Jahr darin, aus einem Stereo-Bildpaareine dritte virtuelle Ansicht zu generieren. Dabei soll der Blickwinkel der virtuellen Ansichtzwischen den beiden realen Ansichten liegen und durch einen Prozentwert frei bestimmbarsein. Das Programm soll in Matlab ohne die Hilfe spezialisierter Toolboxen erstellt werden.

_GUI_:

- GUI ist beschränkt auf die Beispiel Bildpaare

1. GUI aufrufen durch: challengeGUI.m
2. Häckchen für die gewünschten Bildpaare setzte (L1-R1 und L2-R2) 
3. Auf generieren klicken

_OHNE GUI_: 

1. Zeile 18+19: Bildpaare können in challenge.m geändert werden - L1-R1 oder L2-R2
2. Zeile 21: p Faktor ändern
2. Zeile 23: Ki=K2\_opt bzw. K1\_opt für Bildpaar2 bzw 1
3. Zeile 24: Depth Map Bildpaar1: DepthMap=1 
             Depth Map Bildpaar2: DepthMap=2

_EIGENES BILDPAAR OHNE GUI_:

1. Zeile 18+19: Eigene Bildpaare Laden
2. Zeile 21: p Faktor ändern
2. Zeile 23: Eigene Kalibrierungsmatrix 
3. Zeile 24: DepthMap=3



_DISCLAIMER DEPTH MAP_:

Die Depth Map wird für DepthMap=1 und DepthMap=2 nicht live berechnet sondern für jedes Bildpaar geladen. Für DepthMap=3 wird sie live berechnet. 

_GEFORDERTE ERGEBNISBILDER FÜR P=1,0.7,0.45,0.2_

In Ordner Ergebnisse!
