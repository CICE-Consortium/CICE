#!/bin/csh -f

if (-e poll_queue.env) then
  source poll_queue.env
endif

# Parse the job IDs from suite.jobs.  This should work for PBS, Slurm, or IBM LFS but needs
# to be thoroughly tested (so far only tested on PBS)

# Wait for all jobs to finish
foreach line ("`cat suite.jobs`")
  set job = `echo "$line" | sed  's|^[^0-9]*\([0-9]*\).*$|\1|g'`
  set qstatjob = 1
  if (${job} =~ [0-9]*) then
    while ($qstatjob)
      ${ICE_MACHINE_QSTAT} $job >&/dev/null
      set qstatus = $status
#      echo $job $qstatus
      if ($qstatus != 0) then
        echo "Job $job completed"
        set qstatjob = 0
      else
        echo "Waiting for $job to complete"
        sleep 60   # Sleep for 1 minute, so as not to overwhelm the queue manager
      endif
#      echo $qstatjob
    end
  endif
end

