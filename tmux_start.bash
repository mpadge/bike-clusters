#!/bin/sh
cd ./src/
SESSION="Clusters"

tmux -2 new-session -d -s $SESSION

tmux new-window -t $SESSION:1 -k -n main
tmux send-keys -t $SESSION:1 'vim Utils.h' C-m
tmux send-keys -t $SESSION:1 ':' 'tabe Utils.c++' C-m
tmux send-keys -t $SESSION:1 ':' 'tabe Structures.h' C-m

tmux split-window -h
tmux send-keys -t $SESSION:1 'vim Clusters.h' C-m
tmux send-keys -t $SESSION:1 ':' 'tabe Clusters.c++' C-m
tmux select-pane -t 0

tmux new-window -t $SESSION:2 -n routines1
tmux select-window -t $SESSION:2
tmux send-keys -t $SESSION:2 'vim InOut.h' C-m
tmux send-keys -t $SESSION:2 ':' 'tabe InOut.c++' C-m

tmux split-window -h
tmux send-keys -t $SESSION:2 'vim DataProcessing.h' C-m
tmux send-keys -t $SESSION:2 ':' 'tabe DataProcessing.c++' C-m
tmux select-pane -t 0

tmux new-window -t $SESSION:3 -n routines2
tmux select-window -t $SESSION:3
tmux send-keys -t $SESSION:3 'vim Calculations.h' C-m
tmux send-keys -t $SESSION:3 ':' 'tabe Calculations.c++' C-m

tmux split-window -h
tmux select-pane -t 0

cd ..
tmux new-window -t $SESSION:4 -n makefile
tmux select-window -t $SESSION:4
tmux send-keys -t $SESSION:4 'vim makefile' C-m

tmux split-window -h

cd ./src/
# The following command does not work with vim-R
#tmux set-window-option -g automatic-rename off
tmux new-window -t $SESSION:5 -n R
tmux select-window -t $SESSION:5
tmux send-keys -t $SESSION:5 'vim cluster_significance.R' C-m
tmux send-keys -t $SESSION:5 ',rf'

tmux attach -t $SESSION
