#--cal_REWARD
# Gives water from right & left with sound manually. ( Press Space)
# LAST EDIT: HA 12/1/2018

#-------SDL HEADER----------------------------------------------------------------------------------------#
active_buttons = 5;	# Needed for other .SCE files in the experiment (see Settings>Response>Active Buttons).
no_logfile = true;    # Don't write a log file.
begin;

sound {
	wavefile { filename ="SAM_5k_20FM_AM50.wav"; preload = true; };
} soundCue;

# SDL EVENTS
trial{
trial_type = correct_response;
trial_duration = forever;
nothing {};
target_button = 4,5; # Spacebar
}waitForUser;

trial {
	#all_responses = false;	#first_response, but ignore the responses before stimulus_time_in
	trial_type = first_response;
	trial_duration = 2000;	
   sound soundCue;	# Go Cue
	target_button = 2,3;   
} responseWindow;

#-------PCL-----------------------------------------------------------------------------------------------#
begin_pcl;

# SET TASK PARAMETERS AND INITIALIZE TRIAL VARIABLES
include "PR_Parameters.pcl";
int txt =0;
# SETUP OUTPUT PORT
output_port port = output_port_manager.get_port(1);

# PARAMETER WINDOW
parameter_window.remove_all();
int WaterIdx = parameter_window.add_parameter("Water Deliveries");
display_window.draw_text("Press SPACE to test valves, ESC to quit.");
# BEGIN
loop int i until i!=i # Boolean dummy expression
begin	
	waitForUser.present(); # Wait for spacebar press
	txt = response_manager.last_response();
	display_window.draw_text(string(txt) + "Last: " + string(response_manager.last_response()));
	if txt == 4 then;	
		port.send_code(portcode_left, waterAmount_left); #LEFT VALVE
		responseWindow.present();
		stimulus_data last = stimulus_manager.last_stimulus_data();
		display_window.erase();
		display_window.draw_text("Left Given");
		if response_manager.response_count()>0 then
      display_window.draw_text("Left Given Response:" + string(response_manager.last_response()) + " RT:" + string(last.reaction_time()));
		end;
      
   elseif txt==5 then;
      port.send_code(portcode_right, waterAmount_right); #RIGHT VALVE
		parameter_window.set_parameter(WaterIdx,string(i) );
		responseWindow.present();
		stimulus_data last = stimulus_manager.last_stimulus_data();
		display_window.erase();
		display_window.draw_text("Right Given");
		if response_manager.response_count()>0 then
      display_window.draw_text("Right Given Response:" + string(response_manager.last_response())+ " RT:" + string(last.reaction_time()));
		end;
   end;
		i = i+ 1;  
end; # End main loop