############################################################################################
# module to encode in BETA information contained in the XYZ file as additional fields      #
# if molecule I has been loaded from file.xyz where the Jth column contains the property   #
# to plot, one should call beta_load J I to scan the file and store the data into a work   #
# array.                                                                                   #
# from there on, every time a frame is selected, the BETA value for the atoms will be set  #
# to the one read from the file. update can be forced with beta_set_all.                   #
############################################################################################
array set xprop {}
array set xnfr {}
array set xnat {}

#set beta value for a molecule (DOES NOT CHECK IF DATA ACTUALLY EXISTS!)
proc beta_set {args} { 
global xprop
global xnfr
global xnat
lassign $args molid
if { $molid == {} } { set molid [molinfo top] }

set sel [atomselect $molid  "all"] 
set frame [molinfo $molid get frame]
set nat [molinfo $molid get numatoms]
if { $frame >= $xnfr($molid) } { 
  puts "Extended data not loaded for frame $frame"
  return
}
if { $nat != $xnat($molid)} {
  puts "Atoms number in molecule and in extended properties mismatch."
  return
}

puts "setting frame data for mol. $molid"
$sel set beta $xprop($molid,$frame)
$sel delete
}

#sets beta value for all molecules which have been loadad
proc beta_set_all {args} { 
global xnfr
set lmol [molinfo list]
foreach molid $lmol {
  if { [ info exists xnfr($molid) ] > 0 } { beta_set $molid }
}
}


#loads the selected field for the given molecule id (default 4 and "top")
proc beta_load {args} {
global xprop
global xnfr
global xnat
lassign $args field molid
if { $field == {} } { set field 4 }
if { $molid == {} } { set molid [molinfo top] }

#gets molecule info
set fname [molinfo $molid get filename]
set nframes [molinfo $molid get numframes]
set nat [molinfo $molid get numatoms]
set xnfr($molid) $nframes
set xnat($molid) $nat

puts "* reading extra field $field from $fname (molecule $molid)"
set myfile  [ open $fname ] 

for {set iframe 0}  {$iframe< $nframes} { incr iframe } {
  #reads frame by frame and prepares an array with the 
  #extra property (col. 5 in XYZ file) for each frame
  set xprop($molid,$iframe) {}	
  puts "* reading frame $iframe"
  set line [ gets $myfile ]
  if { $line != $nat } {
    puts "Format error in reading XYZ* file"
    return
  }
  set line [ gets $myfile ]
  for { set at 0 } {$at<$nat} {incr at} {
    set line [ gets $myfile ]
    lappend xprop($molid,$iframe) [ lindex $line $field ]
  }
}
close $myfile
beta_set $molid
}

#automatically update beta on frame change
trace variable vmd_frame(0) w beta_set_all
