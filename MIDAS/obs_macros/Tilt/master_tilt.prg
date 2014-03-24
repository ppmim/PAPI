! combine all tilt measurements into one master-tilt_frame

! if new tilt_frame has an appropriate value, write this value in master_tilt_frame


replace/image master_tilt_frame master_tilt_frame tilt_frame/50000,150000=tilt_frame
