#!/bin/bash

#NOTE: To start the servers, we need previusly to start the Pyro nameserver (pyro-ns) from the command-line
python ./server1.py &
python ./server2.py &
python ./server3.py &
python ./server4.py &

