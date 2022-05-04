#include <iostream>
#include <fstream>
#include "tmfe.h"
#include "midas.h"
#include "mvodb.h"
using namespace std;


//Regular Equipment class for MIDAS
class Equipment : public TMFeEquipment
{
public: 
    Equipment(const char* eqname, const char* eqfilename)
      : TMFeEquipment(eqname, eqfilename)
    {
    }
    
};


//Frontend that may be used to link to MIDAS
/*
class Frontend : public TMFrontend
{
public:
    Frontend(){
        FeSetName("Random_Numbers");
        FeAddEquipment(new Equipment("Test Equipment", __FILE__));

    }
    
    TMFeResult FeInit(const std::vector<std::string>& args){
        return TMFrontend::FeInit(args);
    }
    
    //TMFeResult HandleFrontendInit(const std::vector<std::string>& args){
      //  return TMFeOk();
    //};
    
    //void HandleFrontendExit(){};
   
};
*/



//Main function that contains the pointer for connection to MIDAS, and the addition of the equipment.
int main(int argc, char* argv[]){    
TMFE* tfe = TMFE::Instance(); //Here is the pointer that into TMFE, or MIDAS
TMFeResult result;

//Checking to make sure that the pointer was created correctly
if (tfe == NULL){
    cout << "The pointer isnt working" << endl;
    return 1;
}

//-------------------------------------------------------

//Here the program is connecting to the midas framwork. Better yet, its connecting to the server (mhttpd) and the odbedit. 
result = tfe->Connect("Random_Numbers", NULL, "First_Exp");

//Trouble shooting to make sure that the connection was successful.
if(result.error_flag == true){
	cout << "An Error occured when trying to connect to the experiment (MIDAS)" << endl;
    cout << result.error_message << endl;
    return 2;
}
else if(result.error_flag == false){
	cout << "Success in connecting to MIDAS" << endl;
}

//Here is some code for the frontend to read the ODB, and/or create/open a directory and write something to it. Say a file of some name. Then open the file up.


//This creation of equipment is done independently of the classes defined above.

std::vector<std::string> eq_args;
TMFeEquipment* eq = new Equipment("New Equipment", __FILE__);

eq->fEqConfWriteEventsToOdb = true;

eq->fEqConfBuffer = "System 2.0";
eq->fEqConfEventID = 23;
eq->EqInit(eq_args);
eq->EqSetStatus("Starting...", "white");

tfe->AddRpcHandler(eq);
tfe->RegisterRPCs();

	
//Here a key Period Millisecond is created with a value in the ROOT tree.
int fEqConfPeriodMillisec = 5000;
tfe->fOdbRoot->RI("Period Millisecond", &fEqConfPeriodMillisec, true);




//This one uses the classes from above.
//Frontend frontend;
//frontend.FeMain(argc, argv);



//To make sure that the ODB doesnt kick us out after some amount of time.
tfe->SetWatchdogSec(0);


//Here are some conditions on whether the user wants to exit the ODB.
int condition = 1;
string exit;
while (condition){
        cout << "Do you want to exit the Frontend? (Y/n)" << endl;
        cin >> exit;
    if (exit == "Y"){
        
        result = tfe->Disconnect();
        //cout << typeid(result).name() << endl;

        if(result.error_flag == false){
            cout << "Disconnection from MIDAS Successful." << endl;
            return 3;
            
        }
        else if(result.error_flag == true){
            cout << "Disconnection from MIDAS Failed" << endl;
        }
        condition = 0;
    }

}



return 0;


}

