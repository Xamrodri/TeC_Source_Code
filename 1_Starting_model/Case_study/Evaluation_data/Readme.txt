The Parlung No.4 set-up was done as follow:
- DEM: Used the DEM produced by Tom after his downscaling of ERA5-Land for this site
- Glacier outlines: RGI 6.0 (Ndrive), glaciers within the zone of interest were selected manually in QGIS and saved into a shapefile.
- Glacier thickness: Used from Farinotti et al. 2019 (Ndrive). The individual glacier files were extracted using the script 
"Import_Farinotti_icethick.m" from N:\gebhyd\8_Him\Remote_Sensing_Data\Regional_Global\Farinotti2019_composite_thickness_RGI60-all_regions
and then merged into a composite .tif file in QGIS using the merge function.
- Debris extent: Used from Sherler et al. 2018 (Ndrive). Clipped within the extent of glacier outlines in QGIS.
- Debris thickness: Provided by Michael McCarthy, based on the glacier outlines described above.
- Initial snow depth: to be done
- Land cover class: Clipped from the 20° x 20° PROBAV 100m land cover discrete classification map, which can be downloaded at this link:
https://lcviewer.vito.be/download . The clipping was done in QGIS, with a (generously wide, bounding box encompassing our
catchment of interest).
- Soil properties: Based on SOILGRIDS. The soil depth was extracted from the product: depth to bedrock (R horizon) up to 200 cm predicted 
using the global compilation of soil ground observations , which can be downloaded at this link:
https://data.isric.org/geonetwork/srv/eng/catalog.search#/metadata/bfb01655-db81-4571-b6eb-3caae86c037a
The clay, sand and soil organic content per soil layers from SOILGRIDS were downloaded using the R Script written by Pascal Buri
"download_SoilGrids_v2.R". Only the bouding box (lat/lon) and domain name needed to be adapted.
The soil layers were then merged using the R Script from Pascal: "derive_avg_SOMA-SOILGRIDS_soiltexture_maps_v2", which I had to 
adapt slightly (file paths, file names). Note that the script needs to be run one time for each soil texture.
- The catchment outlet was selected based on the location of the stream gauging station. 