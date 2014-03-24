import logging

logging.basicConfig()

log = logging.getLogger("MyApp")
log.setLevel(logging.DEBUG) #set verbosity to show all messages of severity >= DEBUG
log.info("Starting my app %s %d %d", "Number", 1, 2)


