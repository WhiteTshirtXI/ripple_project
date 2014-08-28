function tview(img, range)
%Return a modified gray scale that is similar to the one used in the tview
%program.
%
%Example:
%>> gmap = create_new_gray_scale();
%>> imagesc(image, [0 1000]);
%>> colormap(gmap);
red = cat(3, ones(size(img)), zeros(size(img)), zeros(size(img)));
green = cat(3, zeros(size(img)), ones(size(img)), zeros(size(img)));
r = imagesc(red); hold on
g = imagesc(green);
h = imagesc(img, [range(1), range(2)]); hold off
%colormap(gray)

set(r, 'AlphaData', (img <= 0))
set(g, 'AlphaData', (img >= 65000))

% Uncomment this line to display both red and green pixels
%set(h, 'AlphaData', opacity*(img > 0 & img < 65000))

% Uncomment this line to display only green pixels
set(h, 'AlphaData', (img < 65000))
