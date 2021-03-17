## get_emp_color() retrieves the standardized color for each Waimea metadata category
# If cat_name is specified, returns a single value. If return_all is set to TRUE, returns a vector of coors. 
# arguments:
# cat_name = the name of the category for which a color should be returned (e.g. 'Free-living')
# return_all = if set to TRUE, returns a full named vector of all colors for specified category type(s)
# cat_type = the type of category, either "emp", "habitat", or trophic", useful when combined with return_all

get_waimea_colors <- function(cat_type = c("emp","habitat","trophic","site")){

  # colors as defined in Thompson et al 2017
  # https://github.com/biocore/emp/blob/master/code/colors-styles/empcolors.py
    
  emp_colors <-c("Other" ='#929591', # grey for null value
    'EMP sample'= '#929591', # 'grey'
    'Host-associated'= '#fb9a99',
    'Free-living'= '#e31a1c',
    'Animal'= '#b2df8a',
    'Plant'= '#33a02c',
    'Fungus' = "#bee5bf",
    'Non-saline'= '#a6cee3',
    'Saline'= '#1f78b4',
    'Aerosol (non-saline)'= '#d3d3d3', # 'lightgrey'
    'Animal corpus'= '#ffff00', # 'yellow'
    'Animal distal gut'= '#8b4513', # 'saddlebrown'
    'Animal proximal gut'= '#d2b48c', # 'tan'
    'Animal secretion'= '#f4a460', # 'sandybrown'
    'Animal surface'= '#b8860b', # 'darkgoldenrod'
    'Hypersaline (saline)'= '#87cefa', # 'lightskyblue'
    'Intertidal (saline)'= '#afeeee', # 'paleturquoise'
    'Mock community'= '#ff00ff', # 'fuchsia'
    'Plant corpus'= '#7cfc00', # 'lawngreen'
    'Plant rhizosphere'= '#006400', # 'darkgreen'
    'Plant surface'= '#00fa9a', # 'mediumspringgreen'
    'Fungus corpus' = "#bee5bf",
    'Sediment (non-saline)'= '#ffa07a', # 'lightsalmon'
    'Sediment (saline)'= '#ff6347', # 'tomato'
    'Soil (non-saline)'= '#ff0000', # 'red'
    'Sterile water blank'= '#ee82ee', # 'violet'
    'Surface (non-saline)'= '#000000', # 'black'
    'Surface (saline)'= '#696969', # 'dimgrey'
    'Water (non-saline)'= '#000080', # 'navy'
    'Water (saline)'= '#4169e1' # 'royalblue'
    )
  
  # discrete colors for three habitat types
  
  habitat_colors <- c("Riverine" = "#a6cee3", # purple
                      "Marine"   = "#1f78b4", # blue
                      "Terrestrial" = "#fed9a6", # green
                      "Other" = "black" # gray for null value
                      )
  # discrete colors for trophic levels
  # consumers follow a color ramp from yellow to red
  
  trophic_colors <- c("Environmental" = "#fc8d62",
                      "PrimaryProducer" =  "#66c2a5",
                      # consumers
                      "Consumer" = "#beaed4",
                      "Other"= "black" # gray for null value
                      )
  # discrete colors for sites
  # not used frequently, just get the standard rainbow pattern
  site_colors <- rainbow(14)
  names(site_colors) <-c( "Beach",
                          "Confluence",
                          "DrumRoad",
                          "Entrance",
                          "Estuary",
                          "OutsideOcean",
                          "OutsideRocks",
                          "Ridge",
                          "SharksCove",
                          "ThreeTablesBeach",
                          "TwinRocks",
                          "UppersBeach",
                          "WaimeaBay",
                          "Waterfall")
  
  
  # store all colors in a list
  colors_list <- list(emp_colors,
                      habitat_colors,
                      trophic_colors,
                      site_colors)
  
  type_names <- c("emp","habitat","trophic","site")
  names(colors_list) <- type_names
  
  # convert to data.frame
  all_colors <- lapply(
    type_names,
    FUN = function(x) {
      
      type_colors <- colors_list[[x]]
      
      data.frame(
        Name = names(type_colors),
        Value = type_colors,
        Type = rep(x, length(type_colors))
      )
    }
  )
  
  all_colors <- do.call('rbind', all_colors)
  
  
  # return all colors as a named vector
    some_colors <- as.character(
      all_colors[all_colors$Type %in% cat_type, "Value"])
    names(some_colors) <- all_colors[all_colors$Type %in% cat_type, "Name"]

  return(some_colors)
  
}

