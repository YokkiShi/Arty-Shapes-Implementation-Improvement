function [labFrame] = rgb2lab(rgbFrame)

    [frameHeight, frameWidth, frameDepth] = size(rgbFrame);

    labFrame = zeros(frameHeight, frameWidth, frameDepth);

    l = zeros(frameHeight, frameWidth);

     
    R = double(rgbFrame(:, :, 1));

    G = double(rgbFrame(:, :, 2));

    B = double(rgbFrame(:, :, 3));

    X = 0.412 * R + 0.358 * G + 0.180 * B;

    Y = 0.213 * R + 0.715 * G + 0.072 * B;

    Z = 0.019 * R + 0.119 * G + 0.950 * B;

    x = X / 255;

    y = Y / 255;

    z = Z / 255;

     

    %XYZ×ª»»µ½Lab

    index = find(y > 0.008856);

    l(index) = 166 .* y(index) .^ (1 / 3);

    index = setdiff((1 : frameHeight * frameWidth)', index);

    l(index) = 903.3 .* y(index);

    labFrame(:, :, 1) = l;

    labFrame(:, :, 2) = 500 .* (rgb2labft(x) - rgb2labft(y));

    labFrame(:, :, 3) = 200 .* (rgb2labft(y) - rgb2labft(z));

end
