 

#-------SCL HEADER-------#

scenario = "cal_both_valve";
     
active_buttons = 3;	#how many response buttons in scenario
button_codes = 1,2,3;	
target_button_codes = 1,2,3;

#-------SCL-------#

begin;

trial {
	trial_type = fixed;
	trial_duration = 100;	
	nothing {} waterDelivery_event;
	code=8;
	response_active = true; #still record the licks
}waterDelivery;

#-------PCL-------#

begin_pcl;

display_window.draw_text("Calibrating Both LickPort...");

# SETUP OUTPUT PORT
output_port port = output_port_manager.get_port(1);
int portcode_left = 4;
int portcode_right = 8;

# SET PARAMETERS
int pulsewidth_left = 72;  #max seems to be 50 ms, before water runs off lick port
int pulsewidth_right = 74; 
int interpulse_interval = 500;
int numrepeat = 100;

parameter_window.remove_all();
int WaterIndex = parameter_window.add_parameter("Water Deliveries");

loop
	int i = 1
until
	i > numrepeat
begin

	
	port.send_code(portcode_left, pulsewidth_left); #LEFT VALVE
	#port.send_code(portcode_right, pulsewidth_right); #RIGHT VALVE
	
	waterDelivery.present();
	parameter_window.set_parameter(WaterIndex,string(i));
	i = i + 1;
end;

display_window.draw_text("Calibration Complete");
