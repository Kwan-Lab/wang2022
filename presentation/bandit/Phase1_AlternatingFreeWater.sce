##  Free Water  Version 7/27/18
# Alternates from left to right, trainnig stage for banditTask.
# HA, 12/1/18
#-------HEADER PARAMETERS-------#
scenario = "Phase1 - Free Alternating Water";
active_buttons = 3;							#how many response buttons in scenario
button_codes = 1,2,3;	
target_button_codes = 1,2,3;
response_logging = log_all;				#log all trials
response_matching = simple_matching;	#response time match to stimuli
default_all_responses = true;
begin;

sound {
	wavefile { filename ="SAM_5k_20FM_AM50.wav"; preload = true; };
} soundCue;

#-------SDL EVENTS ('TRIALS')-------#
trial {                           #START TRIAL
   trial_type = fixed;
   trial_duration = 100; 
	nothing {};  
	code = 41;
}startTrial;

trial {
	trial_type = fixed;
	trial_duration = 3000;	#3sec to drink LEFT
	nothing {};
	code=5;
} rewardLeft;

trial {
	trial_type = fixed;
	trial_duration = 3000;	#3sec to drink RIGHT
	nothing {};
	code=6;
} rewardRight;

trial {
	trial_type = fixed;
	trial_duration = 3000; 
	nothing {};
	code=7;
} manual;

trial {
	trial_type = fixed;
	trial_duration = 3000;
	nothing {};
	code=75; 
} norewardLeft;

trial {
	trial_type = fixed;
	trial_duration = 3000;
	nothing {};
	code=76; 
} norewardRight;

trial {                      # MISS 
	trial_type = fixed;
	trial_duration = 3000;
	nothing {};
	code=8;
} pause;

trial {                    #Intertrial No lick period 
	trial_type = fixed;
   trial_duration = 3000; 
   nothing{};
	code=90;
} noLicks;

trial {
   save_logfile {
		filename = "temp.log"; 	# use temp.log in default logfile directory
	};									#save logfile during mid-experiment
}quicksave;

#-----
trial {
	#all_responses = false;	#first_response, but ignore the responses before stimulus_time_in
	trial_type = first_response;
	trial_duration = 5000;
	sound soundCue;	# Go Cue
   code = 21; 
   target_button = 2,3;   
} responseWindow;

#--------PCL---------------
begin_pcl;
display_window.draw_text("Phase 1 - Free Water");
term.print("Starting time:");
term.print(date_time());

# PCL subroutines
string scenario =  "phase1_Alternation";
include "setup_Phase3R_Bandit.pcl";       # Setup file (PCL script)
include "sub_arrayMean.pcl";					# For RT
include "sub_rewardDeliveryPR.pcl"			# rewardDelivery(string target_side)
parameter_window.set_parameter(animalIDind,animalID);

#-------------TRIAL STRUCTURE------------------------------
loop
	int i = 0
	
until
	consecMiss >= max_consecMiss
begin
	parameter_window.set_parameter(trialnumIndex,string(i) + " ("+string(block)+")"); # Trial/block in window for current trial

## Start Period - Response Window
	state = "response_window";
	parameter_window.set_parameter(state_Index,state);	
	
	startTrial.present();
	responseWindow.present();

## Response Period - After Response
	side ="none";
	if response_manager.response_count()>0 then   # lick, not miss
		stimulus_data last = stimulus_manager.last_stimulus_data();
		int  RT            = last.reaction_time();
		if 	(response_manager.last_response()==2) then side = "left"; 
				leftRT.add(RT);  int meanRT = arrayMean(leftRT);   
				parameter_window.set_parameter(leftRT_index,string(RT));
				parameter_window.set_parameter(leftMeanRT_index,string(meanRT));
			elseif 	(response_manager.last_response()==3) then side = "right";
				rightRT.add(RT);  int meanRT = arrayMean(rightRT);  
				parameter_window.set_parameter(rightRT_index,string(RT));
				parameter_window.set_parameter(rightMeanRT_index,string(meanRT));
			elseif   (response_manager.last_response()==1) then side = "manual"; 
		end;
		
		if mod(nTrials_hiP_total,2)!=1 && side=="right"  then    #Fixed reward threshold for this stage
		  state="reward";
	      parameter_window.set_parameter(state_Index,state);
			rewardDeliveryPR(side);       #subrountine give water
			nTrials_hiP_total = nTrials_hiP_total+1;
			right_r = right_r +1;
			parameter_window.set_parameter(right_rIndex,string(right_r));
			parameter_window.set_parameter(nTrials_hiP_totalIndex,string(nTrials_hiP_total));
		elseif  mod(nTrials_hiP_total,2)==1 && side=="left"  then    #Fixed reward threshold for this stage
			state="reward";
	      parameter_window.set_parameter(state_Index,state);
			rewardDeliveryPR(side);       #subrountine give water
			nTrials_hiP_total = nTrials_hiP_total+1;
			left_r = left_r +1;	
			parameter_window.set_parameter(left_rIndex,string(left_r));
			parameter_window.set_parameter(nTrials_hiP_totalIndex,string(nTrials_hiP_total));
		elseif side=="manual" then
			state="manual reward";
	      parameter_window.set_parameter(state_Index,state);
			if consecOutcomeL>=consecOutcomeR then
				port.set_pulse_width(waterAmount_right);
				port.send_code(8);		#give water reward to right:8
				manual.present();
				else
				port.set_pulse_width(waterAmount_left);
				port.send_code(4);		#give water reward to left:4
				manual.present();
			end
		else
			state="no reward";
	      parameter_window.set_parameter(state_Index,state);
			if	side=="right" then
				norewardRight.present(); # no reward
				right_no_r = right_no_r + 1;
				parameter_window.set_parameter(right_no_rIndex,string(right_no_r));
			elseif side=="left" then
				norewardLeft.present(); # no reward
				left_no_r = left_no_r +1;
				parameter_window.set_parameter(left_no_rIndex,string(left_no_r));
			end;
	   end;
		consecMiss = 0; 
		parameter_window.set_parameter(consecmissIndex,string(consecMiss));	
	else 
		state="miss pause";
	   parameter_window.set_parameter(state_Index,state);
		pause.present(); #no response --> next trial	
		consecMiss = consecMiss + 1;
		indMiss	  = indMiss + 1;
		parameter_window.set_parameter(consecmissIndex,string(consecMiss));
		parameter_window.set_parameter(indMissIndex,string(indMiss));
	end;
	
## Control Side Lick
	if side=="left" then
   consecOutcomeL = consecOutcomeL +1;
   consecOutcomeR =0; 
   parameter_window.set_parameter(conOutL,string(consecOutcomeL));
   parameter_window.set_parameter(conOutR,string(consecOutcomeR));
   elseif  side=="right" then
	consecOutcomeR = consecOutcomeR +1;
   consecOutcomeL =0;
   parameter_window.set_parameter(conOutL,string(consecOutcomeL));
   parameter_window.set_parameter(conOutR,string(consecOutcomeR));
   end; 

	
## Window updates - details of block
	block_length = block_length+1;
	i       = i+1;	
	n_trial = i;	# total trial number
	parameter_window.set_parameter(block_Index,string(block));           # display high reward side
	parameter_window.set_parameter(switch_count,string(count_switch));   # display # of Switches
	parameter_window.set_parameter(geo_Index,"ii="+string(ii)+" i_geo="+string(i_geo)+"m"+string(m)); # display i_geo
	
	if  left_r+right_r>2 &&  nTrials_hiP_total>2 then
	  parameter_window.set_parameter(hiP_rateIndex,string(nTrials_hiP_total*100/(n_trial-indMiss))+"%");  # HPS preference 
	  parameter_window.set_parameter(re_Index,string(100*(left_r+right_r)/(n_trial-indMiss))+"%");  # ALL Reward Rate over All trials
	end;

	if (i%5) == 0 then		#every 5 trials, save a temp logfile
	  quicksave.present();
	end;

end;

term.print("Finished:");
term.print(date_time());
term.print("\nTotalLeftLick:");
term.print(string(left_no_r+left_r));
term.print("\nTotalRightLick:");
term.print(string(right_no_r+right_r));
term.print("\nWater Amount Left:");
term.print(string(waterAmount_left));
term.print("\nWater Amount Right:");
term.print(string(waterAmount_right));
