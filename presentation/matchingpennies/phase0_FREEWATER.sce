#-- scenario file --#
# when I press spacebar, give it water
# when mouse licks either port, give it water
# count water given by manualfeed or induced by licks
#

scenario = "phase0_bothports";

active_buttons = 3;	#how many response buttons in scenario
button_codes = 1,2,3;	
#target_button_codes = 1,2,3;
# write_codes = true;	#using analog output port to sync with electrophys
response_logging = log_all;	#log all trials
response_matching = simple_matching;	#response time match to stimuli

begin;


trial {
   save_logfile {
		filename = "temp.log"; 	# use temp.log in default logfile directory
	};									#save logfile during mid-experiment
}quicksave;

trial {
   trial_type = fixed;
   trial_duration = 1000;
	nothing {} startexptevent;
   code=5;
}startexpt;

trial {
	trial_type = fixed;
	trial_duration = 500;	#at least 500ms between water
	nothing {} waterrewardexptevent;
	code=3;
	response_active = true; #still record the licks
}waterreward;

trial {
   trial_type = fixed;
   trial_duration = 100; 
	nothing {} interpulseevent; 
}interpulse;

trial {
   trial_type = first_response;
   trial_duration = forever;
	nothing {} waitlickevent;
   code=0; 
}waitlick;

begin_pcl;

term.print("Starting time:");
term.print(date_time());
logfile.add_event_entry(date_time());

display_window.draw_text("Initializing...");

int num_trials = 400;  # user enters initial value in dialog before scenario
preset int waterAmount_left =60;
preset int waterAmount_right =60;
  
;	# #msec to open water valve
int manualfeed=0;
int leftlick=0;
int rightlick=0;

parameter_window.remove_all();
int manualfeedIndex = parameter_window.add_parameter("Manual feed");
int leftlickIndex = parameter_window.add_parameter("Left Lick");
int rightlickIndex = parameter_window.add_parameter("Right Lick");

# set up parallel port for water reward
output_port port = output_port_manager.get_port(1);

display_window.draw_text("Water reward with left lick or right lick or Spacebar...");

loop
	int i = 1
until
	i > num_trials
begin
	startexpt.present();	
	
	waitlick.present();
	if (response_manager.last_response() == 1) then	#if spacebar
		port.set_pulse_width(waterAmount_left);
		port.send_code(4);		#give water reward to left

		port.set_pulse_width(waterAmount_right);
		port.send_code(8);		#give water reward to right

		waterreward.present();
		manualfeed = manualfeed + 1;
		parameter_window.set_parameter(manualfeedIndex, string(manualfeed));
	elseif (response_manager.last_response() == 2) then	#if licking left
		port.set_pulse_width(waterAmount_left);
		port.set_pulse_width(waterAmount_left);
		port.send_code(4);		#give water reward to left
		waterreward.present();
		leftlick = leftlick + 1;
		parameter_window.set_parameter(leftlickIndex,string(leftlick));
	elseif (response_manager.last_response() == 3) then	#if licking right
		port.set_pulse_width(waterAmount_right);
		port.send_code(8);		#give water reward to right
		waterreward.present();
		rightlick = rightlick + 1;
		parameter_window.set_parameter(rightlickIndex,string(rightlick));		
	end;
	i=i+1;
	
	if (i%5) == 0 then		#every 5 trials, save a temp logfile
		quicksave.present();
	end;
end;

startexpt.present();

display_window.draw_text("Free water session has ended.");
term.print("Ending time:");
term.print(date_time());