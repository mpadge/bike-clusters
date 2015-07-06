#!/bin/sh
cd ./src/
SESSION="Clusters"

tmux -2 new-session -d -s $SESSION

tmux new-window -t $SESSION:1 -k -n main
tmux send-keys -t $SESSION:1 'vim ClusterData.h' C-m
tmux send-keys -t $SESSION:1 ':' 'tabe ClusterData.c++' C-m
tmux send-keys -t $SESSION:1 ':' 'tabe ClusterCalculations.h' C-m
tmux send-keys -t $SESSION:1 ':' 'tabe ClusterCalculations.c++' C-m

tmux split-window -h
tmux send-keys -t $SESSION:1 'vim mainNeutral.h' C-m
tmux send-keys -t $SESSION:1 ':' 'tabe mainNeutral.c++' C-m
tmux send-keys -t $SESSION:1 ':' 'tabe mainActual.h' C-m
tmux send-keys -t $SESSION:1 ':' 'tabe mainActual.c++' C-m
tmux select-pane -t 0

tmux new-window -t $SESSION:2 -n routines1
tmux select-window -t $SESSION:2
tmux send-keys -t $SESSION:2 'vim Utils.h' C-m
tmux send-keys -t $SESSION:2 ':' 'tabe Utils.c++' C-m

tmux split-window -h
tmux select-pane -t 0

cd ../
tmux new-window -t $SESSION:3 -n makefile
tmux select-window -t $SESSION:3
tmux send-keys -t $SESSION:3 'vim makefile' C-m
tmux send-keys -t $SESSION:3 ':' 'tabe tmux_start.bash' C-m

tmux split-window -h
tmux select-pane -t 0

cd ./R/
tmux new-window -t $SESSION:4 -n R
tmux select-window -t $SESSION:4
tmux send-keys -t $SESSION:4 'vim cluster_significance.R' C-m
tmux send-keys -t $SESSION:4 ':' 'cd ../sortpart/R/' C-m
tmux send-keys -t $SESSION:4 ':' 'tabe get-members.R' C-m
tmux send-keys -t $SESSION:4 ':' 'tabe get-clusters.R' C-m
tmux send-keys -t $SESSION:4 ':' 'tabe num-clusts.R' C-m
tmux send-keys -t $SESSION:4 ':' 'cd ../../' C-m

tmux select-window -t $SESSION:1

tmux attach -t $SESSION
