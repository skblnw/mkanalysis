#!/bin/bash

# File to track the status of all processed job IDs
status_file="processed_jobs_status.txt"
touch "$status_file"

# Create a directory to store the output files
output_dir="record_jobid"
mkdir -p "$output_dir"

# Use seq to generate a list of job IDs
start=`grep "Running" "$status_file" | awk '{print $2}' | head -n 1`
if [[ -n "$start" ]]; then
    sed -i "/^JobID: $start (Running)/q" "$status_file"
    sed -i "/^JobID: $start /d" "$status_file"
else
    echo "WARNING: No 'Running' job found in the status file. Using the last job ID."
    start=`grep "JobID" "$status_file" | awk '{print $2}' | tail -n 1`
fi
end=`qstat -a | awk '/^[0-9]/{print $1}' | awk -F. '{print $1}' | tail -n 1`
job_ids=($(seq "$start" "$end"))
echo -e "Processing Job IDs: $start - $end ..."

### Helper Functions ###

# Function to format date
format_date() {
  echo $1 | awk '{printf strftime("%b %d, %Y at %H:%M:%S", $1)}'
}

# Function to extract information from the grep output
process_job_info() {
  job_id=$1
  grep_output=$2

  ### Check if the Job ID is Already Processed ###
  if grep -q "^JobID: $job_id " "$status_file"; then
    echo -e "WARNING: the Job ID ($job_id) is Already Processed"
    return
  fi

  ### Check for Incorrect Termination ###
  incorrect_termination=$(echo "$grep_output" | grep -E ";A;")
  if [ -n "$incorrect_termination" ]; then
    echo "JobID: $job_id (Incorrect Termination)" >> "$status_file"
    return
  fi

  # Extract information
  user=$(echo "$grep_output" | grep ";E;" | awk -F"user=" '{print $2}' | awk '{print $1}')
  jobname=$(echo "$grep_output" | grep ";E;" | awk -F"jobname=" '{print $2}' | awk '{print $1}')
  queue=$(echo "$grep_output" | grep ";E;" | awk -F"queue=" '{print $2}' | awk '{print $1}')
  qtime=$(echo "$grep_output" | grep ";E;" | awk -F"qtime=" '{print $2}' | awk '{print $1}')
  start=$(echo "$grep_output" | grep ";E;" | awk -F"start=" '{print $2}' | awk '{print $1}')
  end=$(echo "$grep_output" | grep ";E;" | awk -F"end=" '{print $2}' | awk '{print $1}')
  res_cput=$(echo "$grep_output" | grep ";E;" | awk -F"resources_used.cput=" '{print $2}' | awk '{print $1}')
  res_mem=$(echo "$grep_output" | grep ";E;" | awk -F"resources_used.mem=" '{print $2}' | awk '{print $1}')
  res_vmem=$(echo "$grep_output" | grep ";E;" | awk -F"resources_used.vmem=" '{print $2}' | awk '{print $1}')
  res_walltime=$(echo "$grep_output" | grep ";E;" | awk -F"resources_used.walltime=" '{print $2}' | awk '{print $1}')
  nodename=$(echo "$grep_output" | grep ";E;" | grep -oP "exec_host=\K(\w+/(\d+)(\+)?)+" | awk -F/ '{print $1}')
  cpu_labels=$(echo "$grep_output" | grep ";E;" | grep -oP "exec_host=\K(\w+/(\d+)(\+)?)+" | sed 's/\w\+\// /g' | tr -d '+')
  cpu_count=$(echo "$cpu_labels" | awk '{print NF}')
  gpu_label=$(echo "$grep_output" | grep ";E;" | grep -oP "exec_gpus=\K(\w+-gpu/(\d+)(\+)?)+" | sed 's/\w\+\-\w\+\// /g' | tr -d '+')
  gpu_count=$(echo "$gpu_label" | awk '{print NF}')

  
  ### Check for Correct Termination ###
  job_terminated=$(echo "$grep_output" | grep ";E;")
  if [ -n "$job_terminated" ]; then
    # Create output file for the current job ID
    output_file="$output_dir/job_${job_id}_info.txt"

    # Write formatted information to the output file
    {
      echo "JobID: $job_id"
      echo "User: $user"
      echo "Job Name: $jobname"
      echo "Queue: $queue"
      echo "Resources Used:"
      echo "Node: $nodename"
      echo "CPU Labels: $cpu_labels"
      echo "CPU Count: $cpu_count"
      echo "GPU Labels: $gpu_label"
      echo "GPU Count: $gpu_count"
      echo "Memory: $res_mem"
      echo "Virtual Memory: $res_vmem"
      echo "Hours Used:"
      echo "Queued Time: $(format_date $qtime)"
      echo "Start Time: $(format_date $start)"
      echo "End Time: $(format_date $end)"
      echo "Walltime: $res_walltime"
      echo "CPU time: $res_cput"
    } > "$output_file"
    
    ### Calculate and Add Queued and Run Hours ###
    run_hours=$(echo $res_walltime | awk -F: '{printf "%.2f", $1 + $2 / 60 + $3 / 3600}')
    queue_hours=$(echo "scale=2;($start - $qtime)/3600" | bc -l)
    echo "Run Hours: $run_hours" >> "$output_file"
    echo "Queued Hours: $queue_hours" >> "$output_file"

    # Add the job ID to the status file
    echo "JobID: $job_id (Correct Termination) Submit Time: $qtime" >> "$status_file"
    return
  fi

  ### Check for Running Jobs ###
  job_running=$(echo "$grep_output" | grep ";S;")
  if [ -n "$job_running" ]; then
    echo "JobID: $job_id (Running)" >> "$status_file"
    return
  fi

  ### Check for Unaccounted Jobs ###
  echo "JobID: $job_id (Unaccounted State)" >> "$status_file"
}

### Main Loop to Process Each Job ID ###
for job_id in "${job_ids[@]}"; do
  grep_output=$(grep ";${job_id}.node200" /opt/gridview/pbs/dispatcher/server_priv/accounting/*)
  process_job_info "$job_id" "$grep_output"
done
