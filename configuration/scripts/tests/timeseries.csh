#!/bin/csh

# Check to see if test case directory was passed
if ( $1 == "-h" ) then
  echo "To generate timeseries plots, this script can be passed a directory"
  echo "containing a logs/ subdirectory, or it can be run in the directory with"
  echo "the log files, without being passed a directory."
  echo "Example: ./timeseries.csh ./annual_gx3_conrad_4x1.t00"
  echo "Example: ./timeseries.csh"
  echo "It will pull the diagnostic data from the most recently modified log file."
  exit -1
endif
set basename = `echo $1 | sed -e 's#/$##' | sed -e 's/^\.\///'`

# Set x-axis limits
  # Manuallyl set x-axis limits
#set xrange = 'set xrange ["19980101":"19981231"]'
  # Let gnuplot determine x-alis limits
set xrange = ''

# Determine if BASELINE dataset exists
if ( $1 == "" ) then
set basefile_dir = "IGNORE"
else
source $1/cice.settings
set basefile_dir = "$ICE_BASELINE/$ICE_BASECOM/$ICE_TESTNAME"
endif

if ( -d $basefile_dir ) then
  set num_basefile = `ls $basefile_dir | grep cice.runlog | wc -l`
  if ( $num_basefile > 0 ) then
    set baseline_exists = 1
    foreach file ($basefile_dir/cice.runlog.*)
      set base_logfile = $file
    end
  else
    set baseline_exists = 0
  endif
else
  set baseline_exists = 0
endif

set fieldlist=("total ice area  (km^2)" \
               "total ice extent(km^2)" \
               "total ice volume (m^3)" \
               "total snw volume (m^3)" \
               "rms ice speed    (m/s)" )

# Get the filename for the latest log
if ( $1 == "" ) then
foreach file (./cice.runlog.*)
  set logfile = $file
end
else
foreach file ($1/logs/cice.runlog.*)
  set logfile = $file
end
endif

# Loop through each field and create the plot
foreach field ($fieldlist:q)
  set fieldname = `echo "$field" | sed -e 's/([^()]*)//g'`
  set search = "'$fieldname'\|istep1"
  rm -f data.txt
  foreach line ("`egrep $search $logfile`")
    if ("$line" =~ *"istep1"*) then
      set argv = ( $line )
      set date = $4
      @ hour = ( $6 / 3600 )
    else
      set data1 = `echo $line | rev | cut -d ' ' -f2 | rev`
      set data2 = `echo $line | rev | cut -d ' ' -f1 | rev`
      echo "$date-$hour,$data1,$data2" >> data.txt
    endif
  end
  set format = "%Y%m%d-%H"

  set output = `echo $fieldname | sed 's/ /_/g'`
  set output = "${output}_${basename}.png"

  echo "Plotting data for '$fieldname' and saving to $output"

# Call the plotting routine, which uses the data in the data.txt file
gnuplot << EOF > $output
# Plot style
set style data points

set datafile separator ","

# Term type and background color, canvas size
set terminal png size 1920,960

# x-axis 
set xdata time
set timefmt "$format"
set format x "%Y/%m/%d"

# Axis tick marks
set xtics rotate

set title "$field (Diagnostic Output)" 
set ylabel "$field" 
set xlabel "Simulation Day" 

set key left top

# Set x-axlis limits
$xrange

if ( $baseline_exists == 1 ) \
  plot "data_baseline.txt" using (timecolumn(1)):2 with lines  lw 2 lt 2 lc 2 title \
       "Arctic - Baseline", \
       "" using (timecolumn(1)):3 with lines lw 2 lt 2 lc 5 title "Antarctic - Baseline", \
       "data.txt" using (timecolumn(1)):2 with lines lw 2 lt 1 lc 1 title "Arctic", \
       "" using (timecolumn(1)):3 with lines lw 2 lt 1 lc 3 title "Antarctic"; \
else \
  plot "data.txt" using (timecolumn(1)):2 with lines lw 2 lt 1 lc 1 title "Arctic", \
       "" using (timecolumn(1)):3 with lines lw 2 lt 1 lc 3 title "Antarctic" \
  
EOF

# Delete the data file
rm -f data.txt
if ( $baseline_exists ) then
  rm -f data_baseline.txt
endif
end
