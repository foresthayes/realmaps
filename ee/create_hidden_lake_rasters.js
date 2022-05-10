// ----------------------------------------------------------------------------- 
// Generate water layer for Glacier National Park area
// ----------------------------------------------------------------------------- 


// Specify region of interset
var roi = ee.Geometry.Rectangle([-114.8000, 48.0688, -112.8530, 49.1000]); 

// Specify paramaters for output

// resolution for transform and output
var res = 10

// output projection
var transform = ee.Projection('EPSG:32612');



// Load and filter Sentinel 2 imagery ------------------------------------------ 

// function to mask clouds from S2 data
function maskS2clouds(image) {
  var qa = image.select('QA60');

  // Bits 10 and 11 are clouds and cirrus, respectively.
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;

  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
      .and(qa.bitwiseAnd(cirrusBitMask).eq(0));

  return image.updateMask(mask);
}

// load imagery
var S2col = ee.ImageCollection("COPERNICUS/S2_SR")
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',20))      // reduce clouds
    .map(maskS2clouds)



// Generate layer of water pixels ----------------------------------------------

var water = S2col
  // filter col to the month of August during the year of interest
  .filterDate(ee.Date.fromYMD(2021, 8, 1), ee.Date.fromYMD(2021, 8, 31))
  .map(function(img){return img.select('SCL')})
  // get most common SCL value during time period  
  .mode()
  // Get SCL value equal to 6 (water)
  .eq(6)
  .clip(roi);


// preview on map
Map.addLayer(water, null, 'Water Layer', 1);

// reproject for based on specified transformation and resolution
var water_proj = water.reproject(transform, null, res)

// Save output to drive
Export.image.toDrive({
  image: water_proj,                  // Name of image to export
  description: 's2_water_2021',       // Export name for file
  scale: res,                         // Define scale in meters
  region: roi,                        // Define region to export
  fileFormat: 'GeoTIFF',              // Define desired file format
  maxPixels: 200000000
});



// Generate elevation layer ----------------------------------------------------

var dem = ee.Image("USGS/3DEP/10m")
  .clip(roi).
  // reproject for based on specified transformation and resolution
  reproject(transform, null, res)

// Save output to drive
Export.image.toDrive({
  image: dem,                         // Name of image to export
  description: 's2_dem',              // Export name for file
  scale: res,                         // Define scale in meters
  region: roi,                        // Define region to export
  fileFormat: 'GeoTIFF',              // Define desired file format
  maxPixels: 200000000
});



// Generate median color RGB image --------------------------------------------

var s2_median = S2col
  .median() // Take the median of all the cloud-free pixels to create one image
  .clip(roi);
  
// create natural color image
var s2_color_vis = {bands: ['B2', 'B3', 'B4'], min: 0, max: 3000, gamma: 1.4};

// Add layers to the map
Map.addLayer(s2_median, s2_color_vis, 'S2 Color Image')

// create RGB image for export
var s2_rgb = s2_median
  .visualize({bands: ['B4', 'B3', 'B2'], min: 0, max: 3000})
  .reproject(transform, null, res);

// Save output to drive
Export.image.toDrive({
  image: s2_rgb,                      // Name of image to export
  description: 's2_Aug2021_color',    // Export name for file
  scale: res,                         // Define scale in meters
  region: roi,                        // Define region to export
  fileFormat: 'GeoTIFF',              // Define desired file format
  maxPixels: 200000000
});
