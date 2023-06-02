#!/bin/bash

# Function to convert the job submit time to seconds since Unix Epoch
function convert_time_to_seconds() {
    date -d@"$1" +"%s"
}

# Function to convert seconds since Unix Epoch to days
function convert_seconds_to_days() {
    echo "$1 / (60*60*24)" | bc
}

# Check for the presence of the log file
if [ ! -f "processed_jobs_status.txt" ]; then
    echo "Error: processed_jobs_status.txt not found."
    exit 1
fi

# Check for time period argument
if [ $# -eq 0 ]; then
    echo "Usage: $0 <time_period_in_days>"
    exit 1
else
    time_period="$1"
fi

# Get the current time in seconds since Unix Epoch
current_time=$(date +"%s")

# Initialize earliest job ID to None
earliest_job_id="None"

# Loop through each line of the log file
while read -r line; do
    # Check if the line contains "Submit Time:"
    if [[ $line == *"Submit Time:"* ]]; then
        # Extract the job submit time from the line
        job_submit_time=$(echo "$line" | grep -oP "\d{10}")

        # Convert job submit time to seconds since Unix Epoch
        job_submit_time_seconds=$(convert_time_to_seconds "$job_submit_time")

        # Calculate the time difference between the current time and the job submit time
        time_diff=$(convert_seconds_to_days "$((current_time - job_submit_time_seconds))")

        # Check if the job submit time is within the time period
        if [ "$time_diff" -le "$time_period" ]; then
            # Extract the job ID from the line
            job_id=$(echo "$line" | grep -oP "JobID: \d+")

            # If earliest_job_id is still None or job_id is smaller, update earliest_job_id
            if [ "$earliest_job_id" == "None" ] || [ "${job_id:7}" -lt "${earliest_job_id:7}" ]; then
                earliest_job_id="$job_id"
            fi
            break
        fi
    fi
done < "processed_jobs_status.txt"

# Print the earliest job ID found
if [ "$earliest_job_id" == "None" ]; then
    echo "No jobs found within the last $time_period days."
else
    echo "Earliest job ID found within the last $time_period days: $earliest_job_id"
fi
