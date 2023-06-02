#!/bin/bash

# Generate an array of all dates within the range
dates=()
for i in $(seq 30 -1 0)
do
    dates+=("$(date -d "$i days ago" +%Y%m%d)")
done

# Create an empty string to hold the output from all the grep operations
grep_output=""

# Convert the dates to log file names and grep them if they exist
for date in "${dates[@]}"; do
    log_file="/opt/gridview/pbs/dispatcher/server_priv/accounting/${date}"
    if [[ -f $log_file ]]; then
        grep_output+=$(grep ";E;" "$log_file" | awk -F";" '{print $3}' | awk -F"." '{print $1}')$'\n'
    fi
done

# After all the grep operations, sort and deduplicate the combined output
job_ids=$(echo "$grep_output" | sort | uniq)

# 30 days ago in seconds since epoch for comparison
thirty_days_ago=$(date -d "30 days ago" +%s)

# Define log file for output
log_file="output_log.csv"

# Write headers to the log file
echo "JobID, User, Queue, CPU Count, GPU Count, Run Hours, Queued Hours" > $log_file

qstat -r | grep " R " > list
while read -r line; do
    job_id=$(echo "$line" | awk -F. '{print $1}')
    user=$(echo "$line" | awk '{print $2}')
    queue=$(echo "$line" | awk '{print $3}')

    qstat_output=$(qstat -f "$job_id")
    cpu_count=$(echo "$qstat_output" | grep "Resource_List.neednodes" | awk -F"ppn=" '{print $2}' | awk -F":" '{print $1}')
    gpu_count=$(echo "$qstat_output" | grep "Resource_List.neednodes" | awk -F"gpus=" '{print $2}')
    
    time_field=$(echo "$line" | awk '{print $11}')
    hour=$(echo "$time_field" | awk -F: '{print $1}')
    min=$(echo "$time_field" | awk -F: '{print $2}')
    run_hours=$(bc -l <<< "scale=2; $hour + $min / 60")

    echo "$job_id, $user, $queue, $cpu_count, $gpu_count, $run_hours, 0" >> "$log_file"
done < list

for job_id in $job_ids; do
  job_info_file="record_jobid/job_${job_id}_info.txt"

  if [[ -f $job_info_file ]]; then
    end_time=$(awk -F': ' '/End Time/ {print $2,$3,$4}' $job_info_file | sed 's/at//' | sed 's/,//' )
    start_time=$(awk -F': ' '/Start Time/ {print $2,$3,$4}' $job_info_file | sed 's/at//' | sed 's/,//' )
    user=$(awk -F': ' '/User/ {print $2}' $job_info_file)
    queue=$(awk -F': ' '/Queue:/ {print $2}' $job_info_file)
    cpu_count=$(awk -F': ' '/CPU Count/ {print $2}' $job_info_file)
    gpu_count=$(awk -F': ' '/GPU Count/ {print $2}' $job_info_file)
    run_hours=$(awk -F': ' '/Run Hours/ {print $2}' $job_info_file)
    queued_hours=$(awk -F': ' '/Queued Hours/ {print $2}' $job_info_file)

    end_time_epoch=$(date -d "$end_time" +%s)
    start_time_epoch=$(date -d "$start_time" +%s)

    # Check if end_time is within last 30 days
    if [[ $end_time_epoch -ge $thirty_days_ago ]]; then
      # If start_time is earlier than 30 days ago, adjust run_hours
      if [[ $start_time_epoch -lt $thirty_days_ago ]]; then
        diff_seconds=$((end_time_epoch - thirty_days_ago))
        run_hours=$(echo "scale=2; $diff_seconds / 3600" | bc)
      fi

      # Write to the log file in CSV format
      echo "$job_id, $user, $queue, $cpu_count, $gpu_count, $run_hours, $queued_hours" >> $log_file
    fi
  fi
done
