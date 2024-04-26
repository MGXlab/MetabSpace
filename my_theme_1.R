#Anneli graph

library(tidyverse)
library(extrafont)
#font_import() #this has to be run once and it imports all of the fonts you have in your system directory, it takes a few minutes
#loadfonts(device = "win")

#prepare your theme

#you can specify your fonts parameters and colors
font <- choose_font("Arial")
fontsize <- 20
basecolor <- "black" #"#7F7F7F" #your base color for both font and lines
  
#and you can specify a there that you will apply to all of your graphs
#this allows having the same style for all graphs
#without the need to copy-past and adjust each time
my_theme <-   theme(
  #remove the background of the plot
  plot.background = element_blank(),
  #and from the panel as well
  panel.background = element_blank(),
  #define the width and color of the axis on the plot
  axis.line = element_line(size = 0.5,
                           color = basecolor),
  #if you use plot title you can specify parameters here
  #PS! use plot title only if you send or show the plot on its own
  #for plots on the slide/thesis use slide title and figure caption
  plot.title = element_text(color = basecolor,
                            size = fontsize,
                            face = "bold",
                            hjust = 0.5, vjust = 1),
  #specify the size and style of the text on the plot, e.g. axis title
  text = element_text(family = font,
                      size = fontsize,
                      color = basecolor),
  legend.key = element_blank(),
  strip.background = element_blank(),
  strip.text = element_text(family = font,
                            size = fontsize,
                            color = basecolor),
  #to remove or adjust the position of the legend
  #"none" - is no legend; "top" "bottom", "right", "left";
  #or by coordinates.
  #c(0, 0) corresponds to the "bottom left"
  #and c(1, 1) corresponds to the "top right" position.
  #legend.position = "none",
  #if you have a legend title and text you can specify font size here
  #here it indicates no legend title
  #legend.title = element_blank(),
  legend.text = element_text(family = font,
                             size = fontsize,
                             color = basecolor),
  #specify axis marks text
  axis.text = element_text(family = font,
                           size = fontsize,
                           color = basecolor),
  #remove tick marks
  axis.ticks = element_blank(),
  #define the ratio of x and y axis
  #PS! for scatter plots it needs to be 1!
  #for predicted - measured plots also adjust the ranges!
  aspect.ratio = 1,
  #adjust the position of the axis title
  axis.title.x = element_text(hjust = c(0.5), vjust = c(1)),
  axis.title.y = element_text(hjust = c(0.5), vjust = c(1)),
)

