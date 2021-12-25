Code and data associated with the paper, Wang H, Ortega HK, Atilgan H, Murphy CE, Kwan AC, Pupil correlates of decision variables in mice playing a competitive mixed-strategy game. eNeuro (2022).

———

Requires the following MATLAB toolboxes: Optimization, Communications

- Add all the subfolders in /MatchingPennies and /bandit to the Path in MATLAB

To start portion on matching pennies:
- Change the variable 'root_path' in line 13 of master_MP_pupillometry.m
e.g. /Users/johndoe/Downloads/pupil_MPBandit-main/MP_Bandit_data/MP/pupilData
- Run master_MP_pupillometry.m

To start portion on two-armed bandit:
- Change the variable 'root_path' in line 13 of master_banditpupil.m
e.g. /Users/johndoe/Downloads/pupil_MPBandit-main/MP_Bandit_data/bandit/pupilData
- Change the variable 'root_path' in line 254 of master_banditpupil.m
e.g. /Users/johndoe/Downloads/pupil_MPBandit-main/MP_Bandit_data/bandit/twopupildata
- Run master_banditpupil.m

