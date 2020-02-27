function gray_image = rgb_to_gray(input_image)
    % Diese Funktion soll ein RGB-Bild in ein Graustufenbild umwandeln. Falls
    % das Bild bereits in Graustufen vorliegt, soll es direkt zurueckgegeben werden.

    [x,y,c]=size(input_image);
    % Check column number
    if(c==1)
        % If there is only one column => Grayscale
        gray_image=input_image;
    else
        if (c==3)
        % 3 columns point to a color image.

        % convert to double
        double_input=double(input_image);

        % formula from homework
        gray_image=double_input(:,:,1)*0.299+double_input(:,:,2)*0.587+double_input(:,:,3)*0.114;
        gray_image=uint8(gray_image);
        else
            % somebody is trying to annoy us
            error('The input image has the wrong number of columns and cannot be grayscale (1 column) or a color image (3 columns).');
        end
    end
end
