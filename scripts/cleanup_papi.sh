#!/bin/bash

# Try to remove any PAPI processes. Not always works.
ps aux | grep -ie papi | awk '{print $2}' | xargs kill -9
