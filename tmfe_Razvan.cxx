//
// tmfe_Razvan.cxx
//
// Simple test for a tmfe c++ frontend with a periodic equipment
//

#include <stdio.h>
#include <signal.h> // SIGPIPE
#include <assert.h> // assert()
#include <stdlib.h> // malloc()
#include <math.h> // M_PI

#include "midas.h"
#include "tmfe.h"

class RazvanEquipment : public TMFeEquipment {

public:
    
    RazvanEquipment(const char * eqname, const char * eqfilename) : TMFeEquipment(eqname, eqfilename) {
    } // constructor doing nothing but inheriting from TMFeEquipment
    
    ~RazvanEquipment() {
    } // destructor doing nothing
    
    TMFeResult HandleRpc(const char * cmd, const char * args, std::string & response) {
       fMfe->Msg(MINFO, "HandleRpc", "RPC cmd [%s], args [%s]", cmd, args);
       
       return TMFeOk();
    } // handles RPC request by doing nothing exceot reporting what it was
    
    TMFeResult HandleBeginRun(int run_number) {
        fMfe->Msg(MINFO, "HandleBeginRun", "Begin run %d!", run_number);
        EqSetStatus("Running", "#00FF00");
        
        return TMFeOk();
    } // handles begging of run by simply bla bla bla 
    
    TMFeResult HandleEndRun(int run_number) {
        fMfe->Msg(MINFO, "HandleEndRun", "End run %d!", run_number);
        EqSetStatus("Stopped", "#00FF00");
        
        return TMFeOk();
    } // handles ending of run by simply bla bla bla 
    
    void HandlePeriodic() {
        double t = TMFE::GetTime(); // gets current time
        double data = 100.0 * sin(M_PI * t / 60); // calculate a sinus with current time 
        fOdbEqVariables->WD("data", data); // write the value to the ODB in the Variables folder
        char status_buf[256]; // declares a string
        sprintf(status_buf, "value %.1f", data); // write value and a bit of bla bla bla into this string
        EqSetStatus(status_buf, "#00FF00"); // write to the system buffer the string "value xxxxx"
    } // this is the work function called periodically when the run is on
};

// main function definition
int main(int argc, char* argv[])
{
    setbuf(stdout, NULL); // I have no clue but is related to printing output
    setbuf(stderr, NULL); // I have no clue but is related to printing errors
    signal(SIGPIPE, SIG_IGN); // set up signals 
    std::vector<std::string> eq_args; // declares a string variable
    TMFE * mfe = TMFE::Instance(); // gets one instance of the MIDAS frontend 
    TMFeResult result = mfe->Connect("tmfe_Razvan", NULL, "First_Exp");
    if (result.error_flag) { // if error write to the console
        fprintf(stderr, "Cannot connect to MIDAS, error \"%s\", bye.\n", result.error_message.c_str());
        
        return 1; // terminates the program
    }
    TMFeEquipment* eq = new RazvanEquipment("tmfe_example", __FILE__); // creates an object of RazvanEquipment type, i.e. one instance of my equipment
    eq->fEqConfWriteEventsToOdb = true;
    
    eq->fEqConfPeriodMilliSec  = 1000; // sets timer to 1 second, i.e. the function HandlePeriodic() will be called by the MIDAS system every second
    eq->fEqConfEventID = 1; // I guess sets the first event ID to one
    eq->fEqConfLogHistory = 1; // logging preference set to one but I don't know what it means
    eq->fEqConfBuffer = "SYSTEM"; // system buffer is used for transmitting events, but in my case I remove all event transmission, I only write to ODB and post a message to the MIDAS logging facility
    eq->EqInit(eq_args); // initialize equipment but not clear to me what the system does exactly at this stage
    eq->EqSetStatus("Starting...", "white"); // write a message to the MIDAS logging facility that the frontend is starting up
    mfe->AddRpcHandler(eq); // RPC registration, declaring the RPC handler
    mfe->RegisterRPCs(); // RPC registration itself
    eq->EqSetStatus("Started...", "white"); // post a message that the frontend is ready to go
    mfe->Disconnect(); // disconnect from MIDAS
    
    return 0; // terminate program
}
