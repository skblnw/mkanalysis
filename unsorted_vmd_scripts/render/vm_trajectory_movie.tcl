# Take Picture 1.0
# ----------------

# REQUIREMENTS: VMD Version 1.0 or greater

# DESCRIPTION:
#   Every time take_picture is called, it takes a picture and saves it to
#   a file using a picture number that is maintained within
#   the body of the procedure. A number of parameters can be sent to
#   take_picture to direct its processing.  In particular, the six
#   available options are
#     (1) reset:  All internal take_picture data is reset to
#           default values
#     (2) frame:  Explicitly specifies the number to be included
#           in the numeric placeholder of the filename
#           specification (see 'format' below)
#     (3) format: Specifies the naming convention to follow
#           when writing output files.  This uses
#           C-like notation to denote placeholders for
#           numbers (see usage examples)
#     (4) method: Specify the type of rendering method to
#           be applied to the current scene.
#     (5) exec:   Specifies a Unix command to evaluate
#     (6) modulo: Outputs a file once only every modulo calls
#           to take_picture.  This is useful when writing
#           only every nth frame in an animation, for
#           instance.


# PROCEDURES:
#   take_picture - a generalized procedure for taking screen shots


# EXAMPLE USAGE:
#   *To write all subsequent output files under the names image.00000.rgb,
#    image.00001.rgb, etc.
#     take_picture format image.%05d.rgb 
#   *To write all subsequent output files in the directory /scratch/disk2
#    under the names snap.0000.rgb, snap.0001.rgb, etc.
#     take_picture format /scratch/disk2/snap.%04d.rgb 
#   *To output a file only once every 5th time take_picture is called
#     take_picture modulo 5 
#   *To use Raster3D rendering instead of snapshot and to output using
#    the naming convention animate.0000.r3d, animate.0001.r3d, etc.
#     take_picture format animate.%04d.r3d
#     take_picture method Raster3d
#     take_picture exec {render < %s -sgi %s.rgb} 
#   *To return all options to default values
#     take_picture reset


# DOWNLOAD THE FILE:
#   take_picture.tcl


# AUTHOR:
#   Andrew Dalke (dalke@ks.uiuc.edu)

# Trajectory Movie Maker 1.0
# --------------------------

# REQUIREMENTS: VMD Version 1.0 or greater

# DESCRIPTION:
#   Type 'make_trajectory_movie_files' to generate a sequence of RGB 
#   snapshots of the VMD Display.  One snapshot will be made for
#   each frame of the animation for your molecule.  For more
#   information on customizing this procedure, see the details
#   regarding take_picture, which are discussed here.
#   For instance, by altering parameters in the call to take_picture,
#   you might choose to save just every fifth frame of the animation.
#   Or, you could specify a rendering option other than snapshot. 
#   A large range of possibilities exist.

#   You can use the SGI program 'mediaconvert' to make an MPEG or
#   a Quicktime movie of the RGB images that you generate.


# PROCEDURES:
#   make_trajectory_movie_files - generates by default one RGB image
#      per animation frame, using the renderer 'snapshot.'  To change
#      these default operations, change parameters in the call to 
#      take_picture in the body of the procedure
#   take_picture - a general, multi-purpose routine for taking
#      VMD screen shots


# DOWNLOAD THE FILE:
#   trajectory_movie.tcl


# AUTHOR:
#   Andrew Dalke (dalke@ks.uiuc.edu)

proc take_picture {args} {
  global take_picture

  # when called with no parameter, render the image
  if {$args == {}} {
    set f [format $take_picture(format) $take_picture(frame)]
    # take 1 out of every modulo images
    if { [expr $take_picture(frame) % $take_picture(modulo)] == 0 } {
      render $take_picture(method) $f
      # call any unix command, if specified
      if { $take_picture(exec) != {} } {
        set f [format $take_picture(exec) $f $f $f $f $f $f $f $f $f $f]
        eval "exec $f"
       }
    }
    # increase the count by one
    incr take_picture(frame)
    return
  }
  lassign $args arg1 arg2
  # reset the options to their initial stat
  # (remember to delete the files yourself
  if {$arg1 == "reset"} {
    set take_picture(frame)  0
    set take_picture(format) "./snap.%5d.bmp"
    set take_picture(method) snapshot
    set take_picture(modulo) 1
    set take_picture(exec)    {}
    return
  }
  # set one of the parameters
  if [info exists take_picture($arg1)] {
    if { [llength $args] == 1} {
      return "$arg1 is $take_picture($arg1)"
    }
    set take_picture($arg1) $arg2
    return
  }
  # otherwise, there was an error
  error {take_picture: [ | reset | frame | format  | \
  method  | modulo ]}
}
# to complete the initialization, this must be the first function
# called.  Do so automatically.
take_picture reset




proc make_trajectory_movie_files {} {
	set num [molinfo top get numframes]
	# loop through the frames
	for {set i 0} {$i < $num} {incr i} {
		# go to the given frame
		animate goto $i
                # force display update
                display update 
		# take the picture
		take_picture 
        }
}