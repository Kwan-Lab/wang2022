
#-- scenario file --#  
# test water valve for delivering water rewards

scenario = "drain_left_valve";

default_background_color = 128, 128, 128;
active_buttons = 3;	#how many response buttons in scenario
button_codes = 1,2,3;	
target_button_codes = 1,2,3;

begin;

begin_pcl;

display_window.draw_text("Giving water reward on LEFT & RIGHT...");

output_port port = output_port_manager.get_port(1);

loop
	int i = 1
until
	i > 1
begin
	port.set_pulse_width(10000);
	port.send_code(4); #LEFT VALVE
	port.send_code(8); #RIGHT VALVE
	
	i = i + 1;
end;
