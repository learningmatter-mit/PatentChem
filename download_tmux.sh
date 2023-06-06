#!/bin/bash

# should be run from root directory of PatentChem repo

for ((i=2001; i<=2023; i++))
do
  # Create a new tmux session
  session_name="dl-$i"
  tmux new-session -d -s "$session_name"

  # Send the command to the session
  tmux send-keys -t "$session_name" "conda activate patents; python download.py --years $i --data_dir /data/patents_data/ --no_uncompress" Enter
done
