# ------------------------------------------------------------------------------
# Create topographic map with satellite imagery overlay for Hidden Lake, GNP
# ------------------------------------------------------------------------------


# Define study area ------------------------------------------------------------

# Bounding box for Hidden Lake, Glacier National Park
x_min = -113.8
x_max = -113.7
y_min = 48.65
y_max = 48.7

bbox <-
  sf::st_bbox(
    c('xmin' = x_min, 'xmax' = x_max, 'ymax' = y_max, 'ymin' = y_min),
    crs = sf::st_crs(4326)
  ) |>
  sf::st_as_sfc() |>
  sf::st_transform(32612) |>
  sf::as_Spatial()



# Load local rasters -----------------------------------------------------------
elev <- 
  raster::raster("./data/hidden_lake_dem.tiff") |> 
  raster::crop(bbox)

# Load color imagery from Sentinel 2
s2_rgb <- 
  raster::stack("./data/hidden_lake_rgb.tiff") |>
  raster::crop(bbox)

water <- 
  raster::raster("./data/hidden_lake_water.tiff") |> 
  raster::crop(bbox) |>                               # crop to study area
  raster::reclassify(matrix(c(0, NA), nrow = 1)) |>   # cnvrt 0 values to NA
  stars::st_as_stars() |>                             # convert to stars
  sf::st_as_sf(as_points = FALSE, merge = TRUE) |>    # convert to sf (polygon)
  smoothr::drop_crumbs(units::set_units(1000, m^2)) |># drop poly < unit area 
  smoothr::smooth(method = 'ksmooth', smoothness = 3) # smooth polygon


# Convert rasters to matrix for use with rayshader -----------------------------

elev_mat <- rayshader::raster_to_matrix(elev)

# FUNCTION to wrangle RGB data to work with rayshader
wrangle_rgb <- function(x){
  # The array must be transposed since rasters and arrays are oriented 
  # differently in R. aperm() is sued to perform a multi-dimensional transpose.
  r_mat = rayshader::raster_to_matrix(x[[1]])
  g_mat = rayshader::raster_to_matrix(x[[2]])
  b_mat = rayshader::raster_to_matrix(x[[3]])
  
  # create array
  rgb_array = array(0, dim = c(nrow(r_mat),ncol(r_mat),3))
  
  rgb_array[,,1] = r_mat/255 #Red layer
  rgb_array[,,2] = g_mat/255 #Blue layer
  rgb_array[,,3] = b_mat/255 #Green layer
  
  rgb_array = aperm(rgb_array, c(2,1,3))
  
  # enhance contrast
  scales::rescale(rgb_array,to=c(0, 1))
  
  return(rgb_array)
}

# Prepare rgb image overlay
rgb_img <- wrangle_rgb(s2_rgb)

# Download OSM data ------------------------------------------------------------

osm_bbox <-
  bbox |> 
  sf::st_as_sfc() |> 
  sf::st_transform(crs = 4326) |> 
  sf::st_bbox()

trails <-
  osmdata::opq(
    bbox = osm_bbox,
    timeout = 25
  ) |> 
  osmdata::add_osm_feature(
    key = "highway", 
    value = c("path", "footway", "steps")
  ) |> 
  osmdata::osmdata_sf() |> 
  {\(x) x$osm_lines}() |> 
  sf::st_transform(32612)

waterway <-
  osmdata::opq(
    bbox = osm_bbox,
    timeout = 25
  ) |> 
  osmdata::add_osm_feature(
    key = "waterway"
  ) |> 
  osmdata::osmdata_sf() |> 
  {\(x) x$osm_lines}() |> 
  sf::st_transform(32612)

hwy <-
  osmdata::opq(
    bbox = osm_bbox,
    timeout = 25
  ) |> 
  osmdata::add_osm_feature(
    key = "highway",
    value = c(
      "motorway",
      "trunk",
      "primary",
      "secondary",
      "tertiary"
    )
  ) |> 
  osmdata::osmdata_sf() |> 
  {\(x) x$osm_lines}() |> 
  sf::st_transform(32612)


# ------------------------------------------------------------------------------

# Specify scale of elevation data (and other rasters)
z_sc = 10

# Create basemap
basemap_rgb <-
  elev_mat |> 
  rayshader::height_shade(
    # white/gray palette to overlay sat imagery on
    texture = (grDevices::colorRampPalette(c("gray60", "#FFFFFF")))(256)
  ) |>
  rayshader::add_shadow(rayshader::ray_shade(elev_mat,lambert=FALSE), 0.7) |>
  rayshader::add_overlay(
    rayshader::sphere_shade(
      elev_mat,
      texture = "bw",
      zscale = z_sc,
      colorintensity = 5
    ),
    alphalayer = 0.6
  ) |>
  rayshader::add_overlay(
    rgb_img,
    alphalayer = .8
  )

# define FUNCTION to overlay features onto basemap
add_osm_features <- function(basemap){
  basemap |> 
  rayshader::add_overlay(
    rayshader::generate_line_overlay(
      waterway, 
      extent = raster::extent(elev),
      linewidth = 1, 
      color = "skyblue2", 
      lty = 1,
      heightmap = elev_mat
    ),
    alphalayer = 0.6
  ) |> 
  # add water layer based on sat imagery
  rayshader::add_overlay(
    rayshader::generate_polygon_overlay(
      water,
      extent = raster::extent(elev),
      heightmap = elev_mat,
      palette = "#4e7982", # imhof3
      linecolor = "#4e7982",
      linewidth = 0
    )
  ) |> 
  # add trails to map
  rayshader::add_overlay(
    rayshader::generate_line_overlay(
      trails, 
      extent = raster::extent(elev),
      linewidth = 2, 
      color = "salmon2", 
      lty = 3,
      heightmap = elev_mat
    ),
    alphalayer = 0.9
  ) |>
  # add highways to map
  rayshader::add_overlay(
    rayshader::generate_line_overlay(
      hwy, 
      extent = raster::extent(elev),
      linewidth = 6, 
      color = "black", 
      lty = 1,
      heightmap = elev_mat
    ),
    alphalayer = 0.9
  ) |> 
  # add center stripe to highway
  rayshader::add_overlay(
    rayshader::generate_line_overlay(
      hwy, 
      extent = raster::extent(elev),
      linewidth = 1, 
      color = "white", 
      lty = 2,
      heightmap = elev_mat
    ),
    alphalayer = 1
  ) 
}


# Visualize shaded 3d map ------------------------------------------------------

feature <- add_osm_features(basemap_rgb)

# 2d view of map
rayshader::plot_map(feature)

# 3d view of map
feature |>
  rayshader::plot_3d(
    elev_mat,
    zscale = z_sc,
    fov = 80,
    theta = 340,
    zoom = 0.6,
    phi = 60,
    windowsize = c(1600, 900)
    , soliddepth = 0
    , background = "gray50"
  )


# Add label for Logan Pass Visitor Center
rayshader::render_label(
  elev_mat,
  long =  300001,
  lat = 5397146,
  altitude = 3500,
  zscale = z_sc,
  extent = raster::extent(elev),
  textcolor = "black",
  textsize = 1.3,
  linecolor = "black",
  linewidth = 2,
  text = "Logan Pass Visitor Center",
  relativez = FALSE,
  clear_previous = TRUE
)

# clear rgl window
rgl::clear3d()

# Visualize NON-shaded map (render_hghquality has it's own shading process) ----

unshaded  <-
  elev_mat |> 
  rayshader::height_shade(
    # white/gray palette to overlay sat imagery on
    texture = (grDevices::colorRampPalette(c("gray60", "#FFFFFF")))(256)
  ) |>
  rayshader::add_overlay(
    rgb_img,
    alphalayer = .8
  ) |> 
  add_osm_features()

# 3d view of map
unshaded |>
  rayshader::plot_3d(
    elev_mat,
    zscale = z_sc,
    fov = 80,
    theta = 340,
    zoom = 0.6,
    phi = 60,
    windowsize = c(1600, 900)
    , soliddepth = 0
    , background = "gray50"
  )


# Add compass
rayshader::render_compass(position = "W", compass_radius = 100)

# Add scalebar
rayshader::render_scalebar(
  position = "S",
  y = 220,
  limits = c(0, 5),
  scale_length = (c(0, 0.66) + 0.01),
  segments = 5,
  label_unit = "km",
  text_switch_side = TRUE
  , offset = 63
  , radius = 8
)

# Add label for Logan Pass Visitor Center
rayshader::render_label(
  elev_mat,
  long =  300001,
  lat = 5397146,
  altitude = 3500,
  zscale = z_sc,
  extent = raster::extent(elev),
  textcolor = "black",
  textsize = 1.3,
  linecolor = "black",
  linewidth = 2,
  text = "Logan Pass Visitor Center",
  relativez = FALSE,
  clear_previous = TRUE
)


# Render 10 AM sunshine
rayshader::render_highquality(
  filename = "./out/hidden_lake_Aug1_10am.png", 
  lightintensity = 850,
  lightdirection = 104,
  lightaltitude = 36,
  text_size = 35,
  line_radius = 2,
  text_offset = c(0, 20, 0), 
  samples = 800,
  scale_text_size = 30,
  width = 1600, 
  height = 900,
  title_text = 'Hidden Lake, August 1st 10:00 AM',
  title_color = "black",
  title_size = 30,
  title_font = "mono"
)



# Animate movement of sun for a single day -------------------------------------

# Generate sequence of times for which to get sun position
tm <- seq(
  lubridate::as_datetime('2022-08-01T00:00:00', tz = "MST"), 
  lubridate::as_datetime('2022-08-01T23:59:59', tz = "MST"), 
  by = "min"
) |> 
  lubridate::with_tz(tzone = "UTC")

# Get sun position for tm sequence
# NOTE: Azmuth is provided with 0 deg = south. East is -0to-180  W is 0to180,
# Convert to Deg. and add 180 to get degrees from North.
sun <-
  suncalc::getSunlightPosition(
  # date= "2022-05-01",
  date = tm,
  lon = -113.8,
  lat = 48.65
) |> 
dplyr::transmute(
  dt = date,
  time = lubridate::as_datetime(dt) |> 
         lubridate::with_tz("MST") |> 
         format("%I:%M %p"),
  minute = lubridate::minute(dt),
  alt_deg = altitude * 180 / pi,
  az_deg = (azimuth * 180 / pi) + 180,
  intensity =
    dplyr::case_when(
      alt_deg < 0 ~ 100,
      alt_deg < 1 ~ 150,
      alt_deg < 5 ~ 250,
      alt_deg < 10 ~350,
      alt_deg < 15 ~450,
      alt_deg < 20 ~550,
      alt_deg < 25 ~600,
      alt_deg < 30 ~650,
      alt_deg < 35 ~700,
      alt_deg < 40 ~750,
      alt_deg < 45 ~800,
      alt_deg > 45 ~850
    )
) |>
dplyr::filter(minute%%5 == 0) |>                      # subset to 5 min int
dplyr::filter(alt_deg > -5)                           # filter to near horizon

# Iterate through sun positions for a given day
for(i in 1:nrow(sun)){
  svMisc::progress(i, max.value = nrow(sun))
  rayshader::render_highquality(
    filename = paste0("./out/tmp/frame_", i, ".png"),
    lightintensity = sun$intensity[i],
    lightdirection = sun$az_deg[i],
    lightaltitude = sun$alt_deg[i],
    text_size = 35,
    line_radius = 2,
    text_offset = c(0, 20, 0),
    samples = 200,
    scale_text_size = 30,
    title_text = sun$time[i],
    title_font = "mono",
    width = 1920,
    height = 1080
  )
}


# Convert .png files into .mp4 or .gif -----------------------------------------

# generate list of frames
png_files <- sprintf("./out/tmp/frame_%d.png", 1:nrow(sun))

# Convert to mp4
av::av_encode_video(png_files, './out/hidden_lake_aug1.mp4', framerate = 30)

# Convert to gif
gifski::gifski(
  png_files,
  gif_file = "hidden_lake_aug1.gif",
  width = 1600,
  height = 900,
  delay = 0.03
)
