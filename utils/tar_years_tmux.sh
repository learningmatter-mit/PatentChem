#!/bin/bash

# should be run from within patents data directory (directory containing folder for each year)

for ((i=2001; i<=2023; i++))
do
  # Create a new tmux session
  session_name="tar-$i"
  tmux new-session -d -s "$session_name"

  # Send the command to the session
  tmux send-keys -t "$session_name" "tar -czf $i.tar.gz $i" Enter
done
