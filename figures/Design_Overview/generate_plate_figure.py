# This script generates an SVG file 96-well plate
# Originally written by Jan-Hendrik Prinz and modified slightly here

# This requires svgwrite to be pip installed

from IPython.display import display_pretty, display_html, display_jpeg, display_png, display_json, display_latex, display_svg
from IPython.display import SVG
import os
from numpy import genfromtxt, rint

wellDiameter = 7;
plateSize = (127.66, 85.48);
plateIndent = 2;
wellLineThickness = 0.5;
platePadding = 20;
cutLinePadding = 4;
leftNodgeSize = 6;
gridThickness = 0.25;
numberingSteps = 1;
holderIndent = 0.74;
rightInShift = 4.09;
paperSize = (135, 80);
placements = {1, 2};
wellSeperation = 0.5;
paperGrabOutlet = 16.91;
wellPlatePadding = 0.2;
markingsOffset = {15, 6};

import svgwrite
import math

svg_document = svgwrite.Drawing(filename = "96-well-plate.svg",
                                size = ("1300px", "900px"))



plate_rows = 8
plate_cols = 12

plate_size_x = 1276
plate_size_y = 855

space = 90
radius = 36

first_x = (1300 - plate_size_x) / 2 + (plate_size_x - space * (plate_cols - 1)) / 2 + 10
first_y = (900 - plate_size_y) / 2 + (plate_size_y - space * (plate_rows - 1)) / 2 + 10

svg_document.add(svg_document.rect(
                                             insert = ((1300 - plate_size_x) / 2,(900 - plate_size_y) / 2),
                                             size = (plate_size_x, plate_size_y),
                                             rx = 20,
                                             stroke_width = "3",
                                             stroke = "black",
                                             fill = "rgb(255,255,255)"
                                             ))

svg_document.add(svg_document.rect(
                                             insert = ((1300 - plate_size_x) / 2 + 20,(900 - plate_size_y) / 2 + 20),
                                             size = (plate_size_x - 40, plate_size_y - 40),
                                             stroke_width = "3",
                                             stroke = "black",
                                             fill = "rgb(255,255,255)"
                                             ))
for row in range(plate_rows):
    svg_document.add(svg_document.text(
                                       text = chr(row + 65),      
                                       insert = (-0.75 * space + first_x, row * space + first_y + 9),
                                       text_anchor = 'middle',
                                       font_size = '40',
                                       alignment_baseline = 'middle',
                                       font_family = 'Futura'
                                       ))

for col in range(plate_cols):
    svg_document.add(svg_document.text(
                                       text = str(col + 1),      
                                       insert = (col * space + first_x, -0.65 * space + first_y),
                                       text_anchor = 'middle',
                                       font_size = '40',
                                       alignment_baseline = 'middle',
                                       font_family = 'Futura'
                                       ))

for row in range(plate_rows):
    for col in range(plate_cols):
        svg_document.add(svg_document.circle(
                                             center = (col * space + first_x, row * space + first_y),
                                             r = radius,
                                            stroke_width = "3",
                                             stroke = "black",
                                             fill = "rgb(255,255,255)"                                             
                                             )
                         )
        
svg_document.save()    
