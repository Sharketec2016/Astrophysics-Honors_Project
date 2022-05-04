/********************************************************************\

  Name:         frontend.c
  Created by:   Stefan Ritt

  Contents:     Experiment specific readout code (user part) of
                Midas frontend. This example simulates a "trigger
                event" and a "periodic event" which are filled with
                random data.
 
                The trigger event is filled with two banks (ADC0 and TDC0),
                both with values with a gaussian distribution between
                0 and 4096. About 100 event are produced per second.
 
                The periodic event contains one bank (PRDC) with four
                sine-wave values with a period of one minute. The
                periodic event is produced once per second and can
                be viewed in the history system.

\********************************************************************/

/**************************************************************************\
 
    Manipulated by: Matthew Buchkowski
    
    Changed Contents: The frontend has been configured to trigger an event, 
    poll that event, and check it against a set criteria. The triggered event 
    is the creation of a random number, between 0 and 1. This random number 
    is then pulled and checked to make sure that it is within a threshold. If 
    the random number, x, is (set limit) < x < 1, then it will be added to a file. 
    The file name is Experiment_(run number#). Also, the number will only ever be 
    able to be written to a file when a run has been started. 
    
    Other contents of the file are connecting, adding, and changing the ODB within 
    the frontend. These portions have been commented out, but are still available.

\***************************************************************************/

#undef NDEBUG // midas required assert() to be always enabled

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h> // assert()
#include <string.h>
#include <time.h>


#include "midas.h"
#include "experim.h"

#include "mfe.h"
#include <iostream>

/*-- Globals -------------------------------------------------------*/

/* The frontend name (client name) as seen by other MIDAS clients   */
const char *frontend_name = "Sample Frontend C Framework";
/* The frontend file name, don't change it */
const char *frontend_file_name = __FILE__;

/* frontend_loop is called periodically if this variable is TRUE    */
BOOL frontend_call_loop = FALSE;

/* a frontend status page is displayed with this frequency in ms */
INT display_period = 3000;

/* maximum event size produced by this frontend */
INT max_event_size = 10000;

/* maximum event size for fragmented events (EQ_FRAGMENTED) */
INT max_event_size_frag = 5 * 1024 * 1024;

/* buffer size to hold events */
INT event_buffer_size = 100 * 10000;

//These values are created for adding new keys to a ODB
int status;
INT k;
BOOL flag;
int size;


//Different file pointers needed for opening up the file
FILE *fp; //File Pointer 
FILE *fp2;

//global variables needed for creating the string, numbers, limit for threshold, and a check if run is executing BOOl
char rand_int[30];
double rand_number; 
double limit = 0.20;
BOOL run_executing = false;



/*-- Function declarations -----------------------------------------*/

INT frontend_init(void);
INT frontend_exit(void);
INT begin_of_run(INT run_number, char *error);
INT end_of_run(INT run_number, char *error);
INT pause_run(INT run_number, char *error);
INT resume_run(INT run_number, char *error);
INT frontend_loop(void);

INT read_trigger_event(char *pevent, INT off);
INT read_periodic_event(char *pevent, INT off);

INT poll_event(INT source, INT count, BOOL test);
INT interrupt_configure(INT cmd, INT source, POINTER_T adr);

//INT run_number(void);

/*-- Equipment list ------------------------------------------------*/

BOOL equipment_common_overwrite = TRUE;

EQUIPMENT equipment[] = {

   {"Trigger",               /* equipment name */
      {1, 0,                 /* event ID, trigger mask */
         "SYSTEM",           /* event buffer */
         EQ_POLLED,          /* equipment type */
         0,                  /* event source */
         "MIDAS",            /* format */
         TRUE,               /* enabled */
         RO_RUNNING |        /* read only when running */
         RO_ODB,             /* and update ODB */
         100,                /* poll for 100ms */
         0,                  /* stop run after this event limit */
         0,                  /* number of sub events */
         1,                  /* don't log history */
         "", "", "",},
      read_trigger_event,    /* readout routine */
   },

   {"Periodic",              /* equipment name */
      {2, 0,                 /* event ID, trigger mask */
         "SYSTEM",           /* event buffer */
         EQ_PERIODIC,        /* equipment type */
         0,                  /* event source */
         "MIDAS",            /* format */
         TRUE,               /* enabled */
         RO_RUNNING | RO_TRANSITIONS |   /* read when running and on transitions */
         RO_ODB,             /* and update ODB */
         1000,               /* read every sec */
         0,                  /* stop run after this event limit */
         0,                  /* number of sub events */
         10,                 /* log history every ten seconds*/
         "", "", "",},
      read_periodic_event,   /* readout routine */
   },

   
//Third piece of test equipment, to show how a new equipment is configured
   {"Third_Test_Equipment",
       {3, 0,                 /* event ID, trigger mask */
         "SYSTEM",           /* event buffer */
         EQ_PERIODIC,        /* equipment type */
         0,                  /* event source */
         "MIDAS",            /* format */
         TRUE,               /* enabled */
         RO_RUNNING |   /* read when running and on transitions */
         RO_ODB,             /* and update ODB */
         1000,               /* read every sec */
         0,                  /* stop run after this event limit */
         0,                  /* number of sub events */
         10,                 /* log history every ten seconds*/
         "", "", "",},
         read_periodic_event,
},
{""}
};










/********************************************************************\
              Callback routines for system transitions

  These routines are called whenever a system transition like start/
  stop of a run occurs. The routines are called on the following
  occations:

  frontend_init:  When the frontend program is started. This routine
                  should initialize the hardware.

  frontend_exit:  When the frontend program is shut down. Can be used
                  to releas any locked resources like memory, commu-
                  nications ports etc.

  begin_of_run:   When a new run is started. Clear scalers, open
                  rungates, etc.

  end_of_run:     Called on a request to stop a run. Can send
                  end-of-run event and close run gates.

  pause_run:      When a run is paused. Should disable trigger events.

  resume_run:     When a run is resumed. Should enable trigger events.
\********************************************************************/

/*-- Frontend Init -------------------------------------------------*/

INT frontend_init()
{
/* put any hardware initialization here */

/*
    //Connecting to the Experiment, and ODB to read/write keys/values to the RO_ODB
   status = cm_connect_experiment("", "First_Exp", "ODB", NULL);
   if(status != CM_SUCCESS){
       std::cout << "Failure to connect to the Experiment" << std::endl;
       return 0;
   }
   else{
       std::cout << "Success in connecting into the Experiment" << std::endl;
        return 1;
   }
   
   
   //For adding a new key, first the ODB and experiment must be linked. This is to prevent confusion for the ODB when logging new keys, as this would corrupt a ODB.
   int result = cm_get_experiment_database(&hDB, NULL);
   size = sizeof(flag);
   db_get_value(hDB, 0, "/Logger/Write data", &flag, &size, TID_BOOL, TRUE);
   
   db_create_key(hDB, 2345, "/Test_key2", TID_INT);
   flag++;
   db_set_value(hDB, 2345, "/Test_key2", &flag, size, 1, TID_INT);
   flag++;
   db_get_value(hDB, 100, "/Test_Key", &flag, &size, TID_BOOL, TRUE);
   
   
   status = cm_disconnect_experiment();
  
   */
   
   /* print message and return FE_ERR_HW if frontend should not be started */
   return SUCCESS;
}

/*-- Frontend Exit -------------------------------------------------*/

INT frontend_exit()
{
    std::cout << "The Frontend is Exiting" << std::endl;
   return SUCCESS;
}

/*-- Begin of Run --------------------------------------------------*/

INT begin_of_run(INT run_number, char *error)
{
   /* put here clear scalers etc. */
    std::cout << "Beginning a Run" << std::endl;
    printf("Run Number is: %d\n", run_number);
    
    //Creating the name of the file. Containing the run number
    char file_name[50] = "/home/lubuntu/online/data/Example_";
    char run[10];
    sprintf(run, "%i", run_number);
    
    strcat(file_name, run);
    fp = fopen(file_name, "w");
        
    //frontend_call_loop = true;
    
    
    
    run_executing = true;
   return SUCCESS;
}

/*-- End of Run ----------------------------------------------------*/

INT end_of_run(INT run_number, char *error)
{
    frontend_call_loop = false;
    run_executing = false;
    fclose(fp);
   return SUCCESS;
}

/*-- Pause Run -----------------------------------------------------*/

INT pause_run(INT run_number, char *error)
{
    frontend_call_loop = false;
   run_executing = false;
    return SUCCESS;
}

/*-- Resuem Run ----------------------------------------------------*/

INT resume_run(INT run_number, char *error)
{
    frontend_call_loop = true;
    run_executing = true;
   return SUCCESS;
}

/*-- Frontend Loop -------------------------------------------------*/

INT frontend_loop()
{
   /* if frontend_call_loop is true, this routine gets called when
      the frontend is idle or once between every event */
   
   
  
   
   return SUCCESS;
}

/*------------------------------------------------------------------*/

/********************************************************************\

  Readout routines for different events

\********************************************************************/

/*-- Trigger event routines ----------------------------------------*/

INT poll_event(INT source, INT count, BOOL test)
/* Polling routine for events. Returns TRUE if event
   is available. If test equals TRUE, don't return. The test
   flag is used to time the polling */
{
    
   int i;
   DWORD flag;
   
   //two if statments. first one checkes to make sure that a run has been started, as to not write to a file otherwise. And the second is the limit condition to check if the random number meets the threshold to be added to the file. Then adding to the file. 
   if(run_executing){
       if(rand_number > limit){
        fputs(rand_int, fp);   
       }
   }
    
   
   for (i = 0; i < count; i++) {
      /* poll hardware and set flag to TRUE if new event is available */
      flag = TRUE;

      if (flag)
         if (!test)
            return TRUE;
   }
   

    
    
   return 0;
}

/*-- Interrupt configuration ---------------------------------------*/

INT interrupt_configure(INT cmd, INT source, POINTER_T adr)
{
    //std::cout << "Interrupt Configureation" << std::endl;
   switch (cmd) {
   case CMD_INTERRUPT_ENABLE:
      break;
   case CMD_INTERRUPT_DISABLE:
      break;
   case CMD_INTERRUPT_ATTACH:
      break;
   case CMD_INTERRUPT_DETACH:
      break;
   }
   return SUCCESS;
}

/*-- Event readout -------------------------------------------------*/

INT read_trigger_event(char *pevent, INT off)
{
    
    //std::cout << "Read Trigger Event" << std::endl;
   UINT32 *pdata;

   /* init bank structure */
   bk_init(pevent);

   /* create structured ADC0 bank */
   bk_create(pevent, "ADC0", TID_UINT32, (void **)&pdata);

   /* following code to "simulates" some ADC data */
   for (int i = 0; i < 4; i++)
      *pdata++ = rand()%1024 + rand()%1024 + rand()%1024 + rand()%1024;
        //std::cout << *pdata << " This is the ADC0 value" << std::endl;

   bk_close(pevent, pdata);

   /* create variable length TDC bank */
   bk_create(pevent, "TDC0", TID_UINT32, (void **)&pdata);

   /* following code to "simulates" some TDC data */
   for (int i = 0; i < 4; i++)
      *pdata++ = rand()%1024 + rand()%1024 + rand()%1024 + rand()%1024;
        //std::cout << *pdata << " This is the TDC0 value" << std::endl;
   bk_close(pevent, pdata);

   /* limit event rate to 100 Hz. In a real experiment remove this line */
   ss_sleep(100);
   
   
  
   //random seed needed for the random number generator
    srand(time(0));
    
    //Here, and string and a double are created. The random number will be assigned to the double, and then transitioned into a string, using the string variable. 
   rand_number = (double)rand()/(double)RAND_MAX;
   
    //sprintf is taking the random number and instead of type casting it, assigning it to a string variable already in memory. 
   sprintf(rand_int, "%f \n", rand_number);
   //std::cout << rand_int << " This is the random number" << std::endl;
    
   
   
   
   //Here i am playing around to see if i could add a new value to the MIDAS bank. 
   //more testing to see how the history plots are being created.
   double rand_number_2 = (double)rand();
   /* create variable length TDC bank */
   bk_create(pevent, "Test_key", TID_UINT32, (void **)&pdata);

   /* following code to "simulates" some TDC data */
   *pdata++ = rand_number_2;

   bk_close(pevent, pdata);
   
   
   
   
   
   

   return bk_size(pevent);
}

/*-- Periodic event ------------------------------------------------*/

INT read_periodic_event(char *pevent, INT off)
{
    //std::cout << "Read Periodic Event" << std::endl;
   UINT32 *pdata;

   /* init bank structure */
   bk_init(pevent);

   /* create SCLR bank */
   bk_create(pevent, "PRDC", TID_UINT32, (void **)&pdata);

   /* following code "simulates" some values in sine wave form */
   for (int i = 0; i < 16; i++)
      *pdata++ = 100*sin(M_PI*time(NULL)/60+i/2.0)+100;

   bk_close(pevent, pdata);
   
   std::cout << pdata << " This is the read periodic event pdata" << std::endl;

   return bk_size(pevent);
}
