%%SPATIAL PLOTS

load("Osat_1.mat")
load("OUTPUT_Kyzylsu_SNOWMAP_13.mat")

%Surface map of the catchment
figure(1)
surf(DTM)

%Spatial SWE
figure(2)
SWE_MAT = reshape(Smelt_spatial_daily, length(x_cell), length(y_cell))
surf(SWE_MAT)

%Spatial SWE
figure(3)
SSN_MAT = reshape(SSN_spatial_daily, length(x_cell), length(y_cell))
surf(SSN_MAT)

%Osat
load("Osat_1.mat")
figure(4)
surf(Osat)

%Osat
load("Ohy_1.mat")
figure(5)
surf(Ohy)

%Matrix
figure(6)
[M, c]=contour(DTM, [1000 1500 2000 2500 3000 3500 4000], "ShowText",true, ...
    "LabelFormat","%d m")
c.LineWidth = 2

%Contour and heatmap overlapped
figure(7)
Ohy_flipped = flipud(Ohy)
imagesc(Ohy_flipped);
caxis([min(Ohy_flipped(:)) max(Ohy_flipped(:))]);
colormap(sky);
colorbar;
set(gca, 'XTick', [], 'YTick', []);

hold on;
DTM_flipped = flipud(DTM)
[M2, c2]=contour(DTM_flipped, [1000 1500 2000 2500 3000 3500 4000 5000], "ShowText",true, ...
    "LabelFormat","%d m", 'LineColor', 'k');
hold off;

%Labels
xlabel({'EAST'})
ylabel({'NORTH'})

%xl = xlim; yl = ylim;

%X ticks
xcord_x = 1:length(x_cell);
xcord_y = 240*ones(1, length(x_cell));
xcord_text = strsplit(num2str(round(x_cell,0)));
xc_selected= 1:20:length(x_cell)

xticks(xc_selected);
xticklabels(xcord_text(xc_selected));

%Y ticks
new_y_cell = flip(y_cell)
ycord_x = 1:length(new_y_cell);
ycord_y = 240*ones(1, length(new_y_cell));
ycord_text = strsplit(num2str(round(new_y_cell,0)));
yc_selected= 1:20:length(new_y_cell)

yticks(yc_selected);
yticklabels(ycord_text(yc_selected));
%{
text(xcord_x(xc_selected), xcord_y(xc_selected),xcord_text(xc_selected) , ...
    'Rotation',45, 'VerticalAlignment','bottom', ...
    'HorizontalAlignment','center');
%}