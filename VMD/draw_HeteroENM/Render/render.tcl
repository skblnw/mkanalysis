display resetview
set string "heteroenmtop.bmp"

scale to 0.035

render Tachyon tmp.dat [format {"C:\Program Files (x86)\University of Illinois\VMD\\tachyon_WIN32.exe"  -mediumshade -trans_max_surfaces 1 tmp.dat -res 2048 1024 -format BMP -o C:/Users/ukevi/Dropbox/QWD/BAR-PH/draw_HeteroENM/Render/%s} $string]

rotate x to -90

set string "heteroenmfront.bmp"

render Tachyon tmp.dat [format {"C:\Program Files (x86)\University of Illinois\VMD\\tachyon_WIN32.exe"  -mediumshade -trans_max_surfaces 1 tmp.dat -res 2048 1024 -format BMP -o C:/Users/ukevi/Dropbox/QWD/BAR-PH/draw_HeteroENM/Render/%s} $string]