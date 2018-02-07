#!/bin/csh -f

# Parse the job IDs from suite.log.  This should work for PBS, Slurm, or IBM LFS but needs
# to be thoroughly tested (so far only tested on PBS)

if (-e suite.jobs) rm -f suite.jobs
set job_id = 0
foreach line ( "`cat suite.log`" )
  if ( $job_id == 1 ) then
    set job_id = 0
    if ( "$line" != " " ) then
      # Grep the job number
      echo "$line" | grep -oP "\d+" | sort -n | tail -1 >> suite.jobs
    endif
  else
    if ( "$line" =~ *'COMPILE SUCCESSFUL'* ) then
      set job_id = 1
    endif
  endif
end

# Wait for all jobs to finish
foreach job ("`cat suite.jobs`")
  while (1)
    ${ICE_MACHINE_QSTAT} $job >&/dev/null
    if ($? != 0) then
      echo "Job $job completed"
      break
    endif
    echo "Waiting for $job to complete"
    sleep 60   # Sleep for 1 minute, so as not to overwhelm the queue manager
  end
end

#rm suite.jobs  # Delete the list of job IDs
