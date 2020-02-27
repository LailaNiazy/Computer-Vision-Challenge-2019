function merkmale = harris_detektor(input_image, varargin)
    % In dieser Funktion soll der Harris-Detektor implementiert werden, der
    % Merkmalspunkte aus dem Bild extrahiert
    %% Input parser aus Aufgabe 1.7
        % In dieser Funktion soll der Harris-Detektor implementiert werden, der
    % Merkmalspunkte aus dem Bild extrahiert
    
    %display('Harris: Inputs');
    
    % Input parser
    % Für jedes Element die Standardwerte und Bedingungen definieren
    valid_input_image = @(x) ismatrix(x);
    
    % segment_length (numerisch, ungerade, > 1): steuert die Groesse des Bildsegments (Standardwert: 15)
    def_sl = 15;
    valid_sl = @(x) isnumeric(x) && isscalar(x) && (x > 1) && (mod(x,2) == 1);
    
    % k (numerisch, k <elementof> [0,1]): gewichtet zwischen Ecken- und Kantenprioritaet (In der Literatur wird oftmals k = 0.05 gesetzt)
    def_k = 0.05;
    valid_k = @(x) isnumeric(x) && isscalar(x) && (x >= 0) && (x <= 1);
    
    % tau (numerisch, >0): legt den Schwellenwert zur Detektion einer Ecke fest (Standardwert: 10^6)
    def_tau = 10^6;
    valid_tau = @(x) isnumeric(x) && isscalar(x) && (x>0);
    
    % do_plot (logical): bestimmt, ob das Bild angezeigt wird oder nicht (Standardwert: False)
    def_dp = false;
    valid_dp = @(x) islogical(x) && isscalar(x);
    
    % min_dist (numerisch, >=1): ist der minimale Pixelabstand zweier Merkmale (Standardwert: 20)
    def_min_dist = 20;
    valid_min_dist = @(x) isnumeric(x) && isscalar(x) && (x >=0);
    
    % tile_size (numerisch): definiert die Kachelgroesse, je nach Eingabe entweder die Seitenlaenge
    % fuer eine quadratische Kachel oder ein Vektor mit zwei Eintraegen fuer Hoehe und Breite. (Standardwert: 200)
    def_tile_size = 200;
    % Validationsfunktion am ende   
    
    % N (numerisch, >=1): ist die maximale Anzahl an Merkmalen innerhalb einer Kachel (Standardwert: 5)
    def_n = 5;
    valid_n = @(x) isnumeric(x) && isscalar(x) && (x >=0);
   
    % alles in einen InputParser laden
    p = inputParser;
    
    addRequired(p, 'input_image', valid_input_image);
    addOptional(p, 'segment_length', def_sl,valid_sl);
    addOptional(p, 'k', def_k, valid_k);
    addOptional(p, 'tau', def_tau, valid_tau);
    addOptional(p, 'do_plot', def_dp, valid_dp);
    addOptional(p, 'min_dist', def_min_dist, valid_min_dist);
    addOptional(p, 'tile_size', def_tile_size, @valid_tile_size);
    addOptional(p, 'N', def_n, valid_n);

    parse(p,input_image,varargin{:});
    
    segment_length = p.Results.segment_length;
    k = p.Results.k;
    tau = p.Results.tau;
    do_plot = p.Results.do_plot;
    min_dist = p.Results.min_dist;
    N = p.Results.N;
    tile_size = p.Results.tile_size;
    
    % tile_size in die richtige form bringen
    if(isscalar(tile_size))
        tile_size=[tile_size tile_size];
    elseif([2,1]==size(tile_size))
        tile_size=tile_size';
    end;
    
    %% Vorbereitung zur Feature Detektion
    % Pruefe ob es sich um ein Grauwertbild handelt
    
    %display('Harris: Preparation');

    [~, ~, x] = size(input_image);
    if(x~=1)
        error("Image format has to be NxMx1");
    end
    
    % Approximation des Bildgradienten
    [Ix,Iy] = sobel_xy(double(input_image));
    
    % Gewichtung
    % Einmal über die Segmentlänge einen Vektor erstellen
    % Je näher zur Mitte der Wert ist, desdo höher
    % Das ganze symmetrisch
    w = min(1:segment_length,segment_length:-1:1);
    % den Vektor normieren
    w = w/sum(w);
  
    % Harris Matrix G
    W=w'*w;
    
    G11=conv2(Ix.*Ix,W,'same');
    G12=conv2(w,w,Ix.*Iy,'same');    
    G22=conv2(Iy.*Iy,W,'same');
    
    %% Merkmalsextraktion ueber die Harrismessung
    %display('Harris: Merkmalsextraktion');
    
    [r,c]=size(G11);
    H=zeros(r,c);

    %Harrismessung
    H=G11.*G22 - G12.*G12 - k*(G11+G22).*(G11+G22);
    
    %Ränder entfernen durch nullsetzen
    edge_size = ceil(segment_length/2);
    H_cleaned = zeros(size(H));
    H_cleaned(edge_size:end-edge_size,edge_size:end-edge_size) = H(edge_size:end-edge_size,edge_size:end-edge_size);
    
    %Werte kleiner tau entfernen
    corners = (H_cleaned>tau).*H_cleaned;
   
    %% Merkmalsvorbereitung
    %display('Harris: Preparation Merkmale');
    % Fuegen Sie einen Nullrand der Breite min_dist um die Matrix corners hinzu. 
    % corners=padarray(corners,[min_dist min_dist]);
    [r,c]=size(corners);
    temp=zeros(r+2*min_dist,c+2*min_dist);
    temp(min_dist+1:end-min_dist,min_dist+1:end-min_dist)=corners;
    corners=temp;
    
    % Indizes
    i=find(corners~=0);
    % Werte (natürlich kein Vergleich gegen 0)
    [~,~,v]=find(corners);
    % Erste Spalte ist das Hauptsortierkriterium, hier Werte einfügen
    sorted=sortrows([v i],'descend');
    % Relevante Indizes entnehmen
    sorted_index=sorted(:,2);
        
    %% Akkumulatorfeld
    %display('Harris: Akkumulation');

    % Das AKKA braucht so viele Felder wie Kacheln ins Bild passen
    [r,c] = size(G11);
    y_tiles = ceil(r/tile_size(1));
    x_tiles = ceil(c/tile_size(2));
    
    AKKA = zeros(y_tiles,x_tiles);
    % merkmale kann muss soviele Merkmale speichern wie vorhanden:
    % 1) Anzahl Kacheln * Anzahl Merkmale pro Kachel
    % 2) Anzahl vorverarbeiteter Merkmale
    [rsi,csi]=size(sorted_index);

    %% Merkmalsbestimmung mit Mindestabstand und Maximalzahl pro Kachel
    %display('Harris: Merkmalsbestimmung');

    % für Übersicht der Merkmalradien 
    unblocked=true(size(corners));
    
    merkmale = [];
    
    c=cake(min_dist);
    
    % über relevante Merkmale iterieren
    for i=sorted_index'
        % checken wo wir sind
        [y,x] = ind2sub(size(corners),i);
        % checken ob das Merkmal schon durch ein anderes geblockt ist
        isfree = unblocked(y,x);
        
        if(isfree)
            % tile finden
            % AKKA deckt nur das ursprüngliche Bild ab
            % Koordinaten im input_image (ohne Padding)
            yii=y-min_dist;
            xii=x-min_dist;
            
            y_tile = ceil(yii/tile_size(1));
            x_tile = ceil(xii/tile_size(2));
            % checken ob noch Platz in der Kachel ist
            if(AKKA(y_tile, x_tile)<N)
                % Merkmal speichern
                merkmale=[merkmale,[xii;yii]];
                
                % Kachelzähler erhöhen
                AKKA(y_tile,x_tile)=AKKA(y_tile,x_tile)+1;
                
                % Merkmalradius eintragen
                % Behandlung der Ränder als Sonderfälle unnötig weil corners bereits um min_dist gepadded ist
                % Grenzen für unblocked Matrix
                uxstart=x-min_dist;uystart=y-min_dist;
                uxend=x+min_dist;uyend=y+min_dist;
               
                %unblocked mit cake verunden
                unblocked(uystart:uyend,uxstart:uxend)=unblocked(uystart:uyend,uxstart:uxend)&c;
            end
        end     
    end
    
    % Plot Routine
    if(do_plot)
        figure
        imshow(input_image);
        hold on
        plot(merkmale(1,:),merkmale(2,:), 'r.');
        hold off
    end;
end

% for the input parser
function valid=valid_tile_size(x)
    if isnumeric(x)
        if isscalar(x)
            valid = true;
        elseif ([1,2]==size(x))
            valid = true;
        elseif ([2,1]==size(x))
            valid = true;
        else
            valid=false;
        end
    else
        valid=false;
    end
end

function [Fx, Fy] = sobel_xy(input_image)
    % In dieser Funktion soll das Sobel-Filter implementiert werden, welches
    % ein Graustufenbild einliest und den Bildgradienten in x- sowie in
    % y-Richtung zurueckgibt.

    % Fuer sobel filter
    s = [1 0 -1; 2 0 -2; 1 0 -1];
    
    Fx = conv2(input_image,s,'same');
    Fy = conv2(input_image,s','same');
end

function Cake = cake(min_dist)
    % Die Funktion cake erstellt eine "Kuchenmatrix", die eine kreisfoermige
    % Anordnung von Nullen beinhaltet und den Rest der Matrix mit Einsen
    % auffuellt. Damit koennen, ausgehend vom staerksten Merkmal, andere Punkte
    % unterdrueckt werden, die den Mindestabstand hierzu nicht einhalten. 
    
    max_dist=2*min_dist+1;
    Cake=false(max_dist);
    % da die 0 nicht mit im array ist
    center=[min_dist+1 min_dist+1];
    for i=1:max_dist
        for j=1:max_dist
            if (norm([i,j]-center)>min_dist)
                Cake(i,j)=true;
            end
        end
    end
end