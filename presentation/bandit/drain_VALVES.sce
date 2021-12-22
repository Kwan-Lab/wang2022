#--drain_VALVES: to clean air bubles in the water tubes. 
# LAST EDIT:  HA, 12/1/18
#-------SDL HEADER----------------------------------------------------------------------------------------#
active_buttons = 5;	# Needed for other .SCE files in the experiment (see Settings>Response>Active Buttons).
no_logfile = true;   # Don't write a log file.
begin;

#-------PCL-----------------------------------------------------------------------------------------------#
begin_pcl;

int pulseWidth = 5000; # Hardcoded; user input not necessary.

output_port port = output_port_manager.get_port(1);
int portcode_left = 4;
int portcode_right = 8;

# DRAIN LEFT AND RIGHT VALVES
display_window.draw_text("Opening LEFT & RIGHT valves for "+string(pulseWidth/1000)+" seconds...");

port.set_pulse_width(pulseWidth);
port.send_code(portcode_left); #LEFT VALVE
port.send_code(portcode_right); #RIGHT VALVE
	
