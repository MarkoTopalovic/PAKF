#pragma once

enum communication_phase
{
	PHASE_DONE			  =-1,
	PHASE_MAIN_FULL       = 1,
	PHASE_SUBSTRUCT_FULL  = 2,
	PHASE_SUBSTRUCT_RHS   = 3,
	PHASE_MAIN_RHS        = 4,
	PHASE_MAIN_VALS 	  = 5
};
