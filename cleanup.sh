#!/bin/bash
ps aux | grep -ie papi.py | awk '{print $2}' | xargs kill -9
