//
//  
//
//  Created by Jen Crodelle on 6/10/2019
// will hold a bunch of functions to be called in a main file
//
//

#ifndef trip_iSTDP_funcs_hpp
#define trip_iSTDP_funcs_hpp

#include <stdio.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <ctime>
#include <sstream>
#include <random>
#include <iomanip>
#include <new>
#include <array>

using namespace std;

double t = 0.0; //time
struct Network
{
    double v, avgVolt, LTD, LTD_cort, backgroundF; // voltage value and the last time step's voltage...
    double o1,r1; // post-synaptic tracer
    double o2, o2_past;
    double iSTDP_trace; //tracer for iSTDP
    int type; // 1 - inhib, 0 - exc
    bool Espike, Ispike; //keeps track of whether neuron spike at given timestep
    
    vector <double> tspN; //spike times of this neuron
    array <double, 1000> synapStrength; //strength of synapse from an input synapse to the cortical cell
    vector < double > tsp_input[1000];  // spike times of the input neuron
    vector < int > synapticConn; //will hold list of LGN SYNAPSES that the neuron receives spikes from

    vector <int> GJConn; //hold GJ connection
    
    double g_excite0, g_inhib0; // conductance
    double g_excite1, g_inhib1; // conductance
    int sisterID;   // sister group number
    // for cortical connections:
    vector <int> corticalConn; // hold list of CORTICAL cells that this neuron SENDS spikes to
    array <double, 400> corticalStrength;   // holds the strength with which this neuron SENDS spikes
    
    double countInputSpikes; //count for the number of spikes received by this neuron from all LGN synapses
    double countExternalSpikes;
    double countExcCorticalSpikes, countInhibCorticalSpikes; //count for number of spikes received by this neuron from other cortical cells

};

// parameter values
double oldv, gE0, gE1,gI0, gI1, t_interval;
double a0,a1,b0,b1,k1,k2,vtilda,k1tilda,k2tilda,tspike;

double vT = -45; //threshold voltage
double vR = -60; //rest voltage
double vReset = -60.0; //reset voltage

double vE =  0.0;
double vI = -80.0;
double gL = 1.0;
double overC = 0.05; // 1/20 ms

// params for connectivity
double sigma_E = 11.0; //ms
double sigma_I = 15.0; //ms
double oversigmaE = 1.0/sigma_E;
double oversigmaI = 1.0/sigma_I;

// params for learning
int N_input = 1000; // number of incoming neurons from LGN
// Pfister params - keep constant
double tau_LTD = 33.7;   //ms
double tau_LTP = 16.8;  //ms
double tau_inhib = 20.0; //ms
double trip_tau_LTD = 114.0; // for o2, tauy
double backgroundNu = 0.5;     //background (sp./ms)
double probSynapseConn = 0.25;  // 20% probability to connect synapse to neuron.
double tau_v = 1000.0;
double overTau_v = 1.0/tau_v;

// for creating intervals:
double r = rand();
double u = r/(double)RAND_MAX;
double newXa = (999.0*u) + 1; // initialize this global variable
int numInhib;
// each LGN synapse needs an r1 tracer
vector <double> r1(1000,0.0); //1000 of them and they start at 0.0
vector <double> r2_past(1000,0.0); //holds the last r2 value
vector <double> r2(1000,0.0); // holds on to r2

// files for information to be exported to:
ofstream finalWeights("LGNweights_final.csv"); //final weights from each LGN synapse to the cortical neurons
ofstream W("weightMatrix_synapses.csv"); // file of LGN weights over time
ofstream W_cortical("weightMatrix_cortex.csv"); // file for Cortical weights over time
ofstream spTimes("neuronSpTimes.csv"); //file for spike times of network
ofstream parameterFile("parameters.txt"); // parameters used in these simulations
ofstream electConn_file("electricConnections.csv"); // electric connections

// function prototypes
// set up the neuron structs
void create_neurons(Network* Neuron, int N, int T,double gmax, int N_input, double backgroundF);
// set up neurons and seed LGN input
void create_neurons_seedLGN(Network* Neuron, int N, int T,double gmax, int N_input, double backgroundF);
// connect neurons with GJs
void connectGJs(Network* Neuron, int N,  int numGJs);
// connect only neighboring neuron with GJs
void connectGJs_neighbors(Network* Neuron, int N,  int numGJs);
// connect GJs with sister cells
void connectGJs_sisters(Network* Neuron, int N, double probGJs);
//connect neurons with synapses using a probability (no structure) and set all initial connections to 0
void setUpZeroSynapses(Network* Neuron, int N,  double probCortConnect);
// give synaptic connections a value.
void connectSynapses(Network* Neuron, int N, double strengthConnect[2][2], double gmax_cortical);
// connect neurons synaptically within some radius (perdiodic in up/down and left/right)
void setUpZeroSynapses_radius(Network* Neuron, int N, int radiusConn[2][2]);
//connect neurons in 1d (just on a line within a radius)
void connect_synap_radius1d(Network* Neuron, int N, int radiusConn[2][2], double strengthConnect[2][2]);
// calculate the input from LGN cells and the resulting plasticity (triplet)
void inputFromLGN(Network* Neuron ,int N, int N_input, double t, double dt, double oversigmaE, double tau_LTD, double A_trip_LTD, double gmax, double flagNewXa,double R0);
// calculate the background drive
void poissonBackground(Network* Neuron, int N, double t, double dt, double oversigmaE, double backgroundNu);
// calculate input from other cotical cells and update plasticity (triplet and iSTDP)
void spikeUpdate(Network* Neuron, int N, double t, double dt, double gmax, double A_trip_LTP_cort,double A_trip_LTP, bool cortexLearnFlag, double A_iSTDP, double tau_inhib, double gmax_cortical, double spikeletSize, double targetRate);
//void runRK(Network* Neuron, int N, int T,  double dt, double gmax, double tau_LTD, double tau_LTP, double tau_inhib,double input_A_tripLTP, double A_trip_LTP_cort, double A_iSTDP, double trip_tau_LTD, int N_input,double t_interval, double backgroundNu, double R0, double targetRate,double gmax_cortical_input, double timeForSynapses);

void create_neurons(Network* Neuron, int N, int T,double gmax, int N_input, double backgroundF)
{
    double pr1,pr2, startValue, countInhib;
    countInhib = 0;
    for (int i = 0; i < N; i++)
    {
        Network n;
        // prob for inhib or exc
        pr2 = (double)rand()/RAND_MAX;
        if (pr2 < 0.2) //inhib
        {
            n.type = 1;
            n.backgroundF = 2.0*backgroundF;
            countInhib = countInhib+1;
            for (int L = 0; L<N_input; L++)
            { n.synapStrength[L] = 0.0;}
            n.sisterID = 0;
        }
        else //make exc and make LGN connections
        {
            n.type = 0;
            n.backgroundF = backgroundF;
            // set sister ID
            n.sisterID = rand()%6+1;
       //     cout << "sister ID = " <<n.sisterID<<endl;
            for (int L = 0; L<N_input; L++)
            {
                pr1 = (double)rand()/RAND_MAX;
                if (pr1 < probSynapseConn)
                {
                    //startValue = 0.5*gmax;
                    startValue = (0.3 + 0.2*((double)rand()/RAND_MAX))*gmax; //random start value for synaptic connection
                    n.synapticConn.push_back(L);
                    n.synapStrength[L] = startValue;
                }
                else
                {
                    n.synapStrength[L] = 0.0;
                }
            }
        }
        
        //all voltages start at rest
        n.v = vR;
        n.o1 = 0.0;//initializes at 0
        n.o2 = 0.0;
        n.r1 = 0.0;//initializes at 0
        n.o2_past = 0.0;
        n.iSTDP_trace = 0.0; // tracer for iSTDP
        
        //initialize conductance to 0
        n.g_excite0 = 0.0;
        n.g_excite1 = 0.0;
        n.g_inhib0 = 0.0;
        n.g_inhib1 = 0.0;
        
        n.avgVolt = 0.0; // will hold the average FR of the cell
        n.LTD = 0.0;
        n.LTD_cort = 0.0;
        //initially no spikes
        n.Espike = false;
        n.Ispike = false;
        // initially no incoming spikes
        n.countInputSpikes = 0.0;
        n.countExternalSpikes = 0.0;
        n.countExcCorticalSpikes = 0.0;
        n.countInhibCorticalSpikes = 0.0;
        Neuron[i] = n; //add this neuron to the network
    }
    cout << "there are " << countInhib << "inhib and " << N-countInhib << " exc neurons " << endl;
    numInhib = countInhib;
}

void create_neurons_seedLGN(Network* Neuron, int N, int T,double gmax, int N_input, double backgroundF)
{
    double pr1, pr2, countInhib,d;
    countInhib = 0;
    for (int i = 0; i < N; i++)
    {
        Network n;
        // prob for inhib or exc
        pr2 = (double)rand()/RAND_MAX;
        if (pr2 < 0.25) //inhib
        {
            n.type = 1;
            n.backgroundF = 2.0*backgroundF;
            countInhib = countInhib+1;
            for (int L = 0; L<N_input; L++)
            { n.synapStrength[L] = 0.0;}
        }
        else //make exc and make LGN connections
        {
            n.type = 0;
            n.backgroundF = backgroundF;
            for (int L = 0; L<N_input; L++)
            {
                pr1 = (double)rand()/RAND_MAX;
                if (pr1 < probSynapseConn)
                {
                    d = L/10.0 - i; // radius of connectivity for LGN to cortex
                    // cout << "i = " << i << ", d = " << d << endl;
                    if (d> 50)
                    {d = 100 -d;}
                    if(d< -50)
                    {d = d+100;}
                    n.synapticConn.push_back(L);
                    n.synapStrength[L] = 0.5*gmax*exp(-0.5*(d/50.0)*(d/50.0));// initial seeding of map
                }
                else
                {
                    n.synapStrength[L] = 0.0;
                }
            }
        }
    
        //all voltages start at rest
        n.v = vR;
        n.o1 = 0.0;//initializes at 0
        n.o2 = 0.0;
        n.r1 = 0.0;//initializes at 0
        n.o2_past = 0.0;
        n.iSTDP_trace = 0.0; // tracer for iSTDP
        
        //initialize conductance to 0
        n.g_excite0 = 0.0;
        n.g_excite1 = 0.0;
        n.g_inhib0 = 0.0;
        n.g_inhib1 = 0.0;
        
        n.avgVolt = 0.0; // will hold the average FR of the cell
        n.LTD = 0.0;
        n.LTD_cort = 0.0;
        //initially no spikes
        n.Espike = false;
        n.Ispike = false;
        // initially no incoming spikes
        n.countInputSpikes = 0.0;
        n.countExternalSpikes = 0.0;
        n.countExcCorticalSpikes = 0.0;
        n.countInhibCorticalSpikes = 0.0;
        Neuron[i] = n; //add this neuron to the network
    }
    cout << "there are " << countInhib << "inhib and " << N-countInhib << " exc neurons " << endl;
    numInhib = countInhib;
}

void connectGJs_neighbors(Network* Neuron, int N,  int numGJs)
{
    int foo, realFoo;
    double pr2;
    // calculate some percentage of GJs
    double probConnectGJ = (double) 0.4;
    int countForGJ = 0;
    // loop over all neurons
    for (int i =0; i<N; i++)
    {   // count number of GJ pairs to stop when real numGJs
        // loop over neighboring cells: (just on either side)
        for (int j=-1; j<2; j=j+2)
        {
            pr2 = (double)rand()/RAND_MAX;
            // connect only two excitatory neurons
            if ((Neuron[i].type == 0) && (Neuron[j].type ==0))
            {
                // connect electrically if within prob and not yourself
                if( (pr2 < probConnectGJ) & (countForGJ < numGJs))
                {
                    foo = i + j;
                    if (foo < 0)
                    {realFoo = N + foo;}
                    else if (foo > N)
                    {realFoo = foo - (N+1);}
                    else
                    {realFoo = foo;}
                    if((Neuron[i].GJConn.size() == 0) && (Neuron[realFoo].GJConn.size()==0))
                    {
                        Neuron[i].GJConn.push_back(realFoo); // put this into the vector of GJ connections
                        Neuron[realFoo].GJConn.push_back(i); // put this into the vector of GJ connections
                        countForGJ = countForGJ+2;
                        cout << "Neuron " << i << " GJ with Neuron " << realFoo << " and are both type " << Neuron[i].type << " and " << Neuron[j].type <<endl;
                    }
                }
            }
        }
    }
}

// just connect GJs in pairs
void connectGJs(Network* Neuron, int N, int numGJs)
{
      cout << "connect GJs in pairs" << endl;
    int foo;
    int countForGJ = 0;
    // loop over all neurons
    for (int i =0; i<N; i++)
    {   // if this neuron is exc and isn't already connected to another neuron
        if ((Neuron[i].type == 0) && (Neuron[i].GJConn.size() == 0))
        {
            // choose another neuron randomly
            foo = rand()%N;
            // if that neuron is exc and isn't already connected to another neuron:
            if((Neuron[foo].type ==0) && (Neuron[foo].GJConn.size() == 0))
            {
                if (countForGJ < numGJs)
                {
                    Neuron[i].GJConn.push_back(foo); // put this into the vector of GJ connections
                    Neuron[foo].GJConn.push_back(i); // put this into the vector of GJ connections
                    countForGJ = countForGJ+2;
                //cout << "Neuron " << i << " GJ with Neuron " << realFoo << " and are both type " << Neuron[i].type << " and " << Neuron[foo].type <<endl;
                }
            }
        }
    }
    cout << "there are " << countForGJ << " many neurons connected by GJ " << endl;
}

void connectGJs_sisters(Network* Neuron, int N, double probGJs)
{
    cout << "connect GJs in sisters" << endl;
    double pr;
    int countForGJ = 0;
    // loop over all neurons
    for (int i =0; i<N; i++)
    {   // if this neuron is exc and isn't already connected to another neuron
        if (Neuron[i].type == 0)
        {
            for (int j =0; j<N; j++)
            { if((Neuron[i].sisterID == Neuron[j].sisterID) && (i != j))
                {
                // prob of connecting to sister cell
                pr = (double)rand()/RAND_MAX;
                if (pr < probGJs)
                {   // connect these cells
                    Neuron[i].GJConn.push_back(j); // put this into the vector of GJ connections
                    Neuron[j].GJConn.push_back(i); // put this into the vector of GJ connections
                    countForGJ = countForGJ+1;
                }
                }
            }
        }
    }
    //gets rid of duplicate connections
    for(int j=0; j<N; j++)
    {
            sort(Neuron[j].GJConn.begin(), Neuron[j].GJConn.end());
            Neuron[j].GJConn.erase(unique(Neuron[j].GJConn.begin(), Neuron[j].GJConn.end()), Neuron[j].GJConn.end());
    }
    cout << "there are around " << countForGJ << " many pairs of GJs " << endl;
}

void setUpZeroSynapses_radius(Network* Neuron, int N, int radiusConn[2][2])
{
    // connects each neuron by +/- radiusSynConn in the up/down and left/right directions
    int X1, Y1, X2, Y2,radiusSynConn;
    double startVal;
    int M = sqrt(N);
    
    for (int i=0; i<N; i++)
    {
        // calculate a new startvalue for EE connections:
        X1 = i%M + 1;
        Y1 = floor(i/M) + 1;
        
        for (int j = 0; j<N; j++)
        {
        if (i != j) // don't connect to yourself and only to those within 40 units
        {   // x,y coordinates of the ith neuron -- the first is (1,1) NOT (0,0)
          //  startVal = 0.7*((double)rand()/RAND_MAX)*gmax_cortical;
          //  strengthConnect[0][0] = startVal;
            // now for jth neuron
            X2 = j%M + 1;
            Y2 = floor(j/M) + 1;
            Neuron[i].corticalStrength[j] = 0.0;
            // find radius based on the cell type (I or E)
            radiusSynConn = radiusConn[Neuron[i].type][Neuron[j].type];
            // if within the radius NOT including periodicity
            if ((abs(X1-X2) < radiusSynConn) && (abs(Y1-Y2) < radiusSynConn))
            {
                Neuron[i].corticalConn.push_back(j); // put this into the vector of cortical connections
              //  Neuron[i].corticalStrength[j] = strengthConnect[Neuron[i].type][Neuron[j].type]; //strength received from neuron j
            }
            else if ((abs(X1-X2) < radiusSynConn)  && (abs(Y1-Y2) > M-radiusSynConn))
            {
                Neuron[i].corticalConn.push_back(j); // put this into the vector of cortical connections
             //   Neuron[i].corticalStrength[j] = strengthConnect[Neuron[i].type][Neuron[j].type]; //strength received from neuron j
            }
            else if ((abs(Y1-Y2) < radiusSynConn) && (abs(X1-X2) > M-radiusSynConn ))
            {
                Neuron[i].corticalConn.push_back(j); // put this into the vector of cortical connections
           //     Neuron[i].corticalStrength[j] = strengthConnect[Neuron[i].type][Neuron[j].type]; //strength received from neuron j
            }
            else if ((abs(X1-X2) > M-radiusSynConn) && (abs(Y1-Y2) > M-radiusSynConn ))
            {
                Neuron[i].corticalConn.push_back(j); // put this into the vector of cortical connections
           //     Neuron[i].corticalStrength[j] = strengthConnect[Neuron[i].type][Neuron[j].type]; //strength received from neuron j
            }
           // else {Neuron[i].corticalStrength[j] = 0.0;}
        }
        else
        {Neuron[i].corticalStrength[j] = 0.0;}
    }
    }
}


//set up synaptic structure, but initialize to 0
void setUpZeroSynapses(Network* Neuron, int N,  double probCortConnect)
{
    //cortical connections start at 0 with all-to-all connectivity
    for (int i=0; i<N; i++)
    {
        for (int j = 0; j<N; j++)
        {
            double pr2 = (double)rand()/RAND_MAX;
            if (i != j) // don't connect to yourself
            {
                if (pr2 < probCortConnect)
                {
                    Neuron[i].corticalConn.push_back(j); // put this into the vector of cortical connections
                }
            }
            // set all connections to 0
            Neuron[i].corticalStrength[j] = 0.0;
        }
    }
}

void connect_synap_radius1d(Network* Neuron, int N, int radiusConn[2][2], double strengthConnect[2][2])
{
    // connects each neuron by +/- radiusSynConn in the up/down and left/right directions
    int radiusSynConn;
    
    for (int i=0; i<N; i++)
    {
        for (int j = 0; j<N; j++)
        {
            if (i != j) // don't connect to yourself and only to those within 40 units
            {
                radiusSynConn = radiusConn[Neuron[i].type][Neuron[j].type];
                if (abs(i-j) < radiusSynConn | abs(i-j) > (N-radiusSynConn))
                {
                    Neuron[i].corticalConn.push_back(j); // put this into the vector of cortical connections
                    Neuron[i].corticalStrength[j]= strengthConnect[Neuron[i].type][Neuron[j].type]; //0.0; //set initial strength to 0 -- strength received from neuron j
                }
                else {Neuron[i].corticalStrength[j] = 0.0;}
            }
            else
            {Neuron[i].corticalStrength[j] = 0.0;}
        }
    }
}

// just put the synaptic strengths in:
void connectSynapses(Network* Neuron, int N, double strengthConnect[2][2], double gmax_cortical)
{
    cout << "start cortical synapses at higher value" << endl;
    int numConn, n1;
    double r1, val;
    //cortical connections start at 0 with all-to-all connectivity
    for (int i=0; i<N; i++)
    {
        r1 = (double)rand()/RAND_MAX;
        numConn = Neuron[i].corticalConn.size();
        for (int j = 0; j<numConn; j++)
        {
            n1 = Neuron[i].corticalConn[j];
            val = strengthConnect[Neuron[i].type][Neuron[n1].type];
           // cout << "type 1 =" << Neuron[i].type << " type 2 = " << Neuron[n1].type << " strength = " << val << endl;
            if ((Neuron[i].type == 0) && (Neuron[n1].type == 0))
            {   //val = gmax_cortical;}
                //val = 0.0;
                val = (0.25 + 0.1*((double)rand()/RAND_MAX))*gmax_cortical; //random start value for synaptic connection
                //cout << val << endl;
                Neuron[i].corticalStrength[n1] = val; // put this into the vector of cortical connections
            }
            else{Neuron[i].corticalStrength[n1] = val;}
        }
    }
}

void inputFromLGN(Network* Neuron ,int N, int N_input, double t, double dt, double oversigmaE, double tau_LTD, double gmax, double flagNewXa,double R0)
{
    double r, u, newu, newr2, tsp, FR, synStrength, updateStrength;
    double R1 = 20.0; //80.0;
    double sigmaRate = 80.0;
    int lengthTsp;
    
    if (flagNewXa == 1) // need to determine a new interval
    {
        r = (double) rand();
        u = r/(double)RAND_MAX;
        ::newXa = (999.0*u) + 1; // number between 0 and 1
        newr2 = (double)rand();
        newu = newr2/(double)RAND_MAX;
        ::t_interval = t-log(newu)*20.0; //calculate the next interval (on avg 20ms in length)
        flagNewXa = 0;
    }
    
    //for each input synapse find spike and update plasticity
    for(int j=0; j<N_input; j++)
    {
        // set spike time to current time (to add later)
        tsp = t;
        while (tsp <= t+dt)
        {
            r = rand();
            u = r/RAND_MAX;
            while (u == 0) //make sure it's not 0, bc then log(0) is inifinite
            {   r = rand();
                u = r/RAND_MAX;}
            // set firing rate of this cell based on its number:
            FR = R0 + R1*(exp(-((::newXa-j)*(::newXa-j))/(2*sigmaRate*sigmaRate)) + exp(-((::newXa+N_input-j)*(::newXa+N_input-j))/(2*sigmaRate*sigmaRate))+exp(-((::newXa-N_input-j)*(::newXa-N_input-j))/(2*sigmaRate*sigmaRate))); // rate determined from Gaussian.
            FR = FR/1000.0; // convert to per ms
            tsp = -log(u)/FR + tsp;
            //if this spike time is still in the right interval of time..
            if (tsp <= t+dt)
            {   // update the tracer for this LGN synapse to say it spiked:
                r1[j] = r1[j] + 1.0;
                r2[j] = r2[j] + 1.0;
                for (int KK = 0; KK<N; KK++)
                {
                    // put this spike time in the list of input spike times for the neurons to which it is connected:
                    // if synapse j is in the list of connections for neuron KK
                    // (these are only exc)
                    if (find(Neuron[KK].synapticConn.begin(), Neuron[KK].synapticConn.end(), j) != Neuron[KK].synapticConn.end())
                    {
                        //put this in the list of spike times for this neuron
                        Neuron[KK].tsp_input[j].push_back(tsp);
                        Neuron[KK].countInputSpikes = Neuron[KK].countInputSpikes + 1;
                        // update the conductance of this neuron for receiving this spike:
                        synStrength = Neuron[KK].synapStrength[j];
                        Neuron[KK].g_excite1 = Neuron[KK].g_excite1 + synStrength*exp(-(t - tsp)*oversigmaE);
                        
                        //if this list gets too big, take one out (for memory issues)
                        if (Neuron[KK].tsp_input[j].size() > 10.0)
                        { Neuron[KK].tsp_input[j].erase(Neuron[KK].tsp_input[j].begin());}
                        
                        //total number of cortical spikes for neuron KK:
                        lengthTsp = Neuron[KK].tspN.size();
                        //if this neuron has spikes, then we take
                        if (lengthTsp > 0)
                        {   //LTD since neuron spiked before LGN
                            updateStrength = -Neuron[KK].o1*Neuron[KK].LTD;
                            Neuron[KK].synapStrength[j] = Neuron[KK].synapStrength[j] + updateStrength*gmax;
                            //cout << "LTD: strength = " << updateStrength<<endl;
                        }
                            // if too small, set to 0
                            if (Neuron[KK].synapStrength[j] < 0.0)
                            {Neuron[KK].synapStrength[j] = 0.0;}
                    } // end if this neuron receives input from this LGN cell
                } // end loop over neurons
            } // end if this LGN spike is within the time step
        } //ends while
    } //ends for
} // ends function

void spikeUpdate(Network* Neuron, int N, double t, double dt, double gmax, double A_trip_LTP_cort, double A_trip_LTP, bool cortexLearnFlag, double A_iSTDP, double tau_inhib, double gmax_cortical, double spikeletSize, double targetRate)
{
    int numGJ, synConN, synConNI, GJcell, foo_synConN, indexL,numLGNconnect,numSynConn, numSynConn_elect;
    double lengthTsp, updateStrength,updateStrength2,updateStrength_iSTDP, Tsp;
    double vT_local = -45.0;
    // loop over all neurons
    for (int l = 0; l<N; l++)
    {
        // if neuron l spiked in this time step
        if (Neuron[l].Espike)
        {
            // update this neuron's tracer:
            Neuron[l].o1 = Neuron[l].o1 + 1.0;
            Neuron[l].o2 = Neuron[l].o2 + 1.0;
            Neuron[l].r1 =  Neuron[l].r1 + 1.0;
            Neuron[l].iSTDP_trace =  Neuron[l].iSTDP_trace + 1.0;
            
            // send spike to post-synaptic cells:
            Tsp = Neuron[l].tspN.back(); //get last spike time for this neuron
            numSynConn = Neuron[l].corticalConn.size(); // get number of synaptic connections
            
            // if the spiking cell is EXC:
            if (Neuron[l].type == 0)
            {   //loop over connections
                Neuron[l].countExcCorticalSpikes = Neuron[l].countExcCorticalSpikes + 1.0; //update spike counter
                for (int ii=0; ii<numSynConn; ii++)
                {
                    // get neuron number of connected cell:
                    synConN = Neuron[l].corticalConn[ii]; //neuron number of iith connection to neuron l
                    Neuron[synConN].g_excite1 = Neuron[synConN].g_excite1 + Neuron[synConN].corticalStrength[l]*exp(-(t - Tsp)*oversigmaE); //update EXC conductance
                    // STDP only if time is right
                    if (cortexLearnFlag)
                    {  //if the connected neuron is exc
                        if (Neuron[synConN].type == 0)
                        {
                            //  plasticity from l --> synConN (ONLY LTD)
                            updateStrength = -Neuron[synConN].o1*Neuron[synConN].LTD_cort; // LTD should be for the neuron who is receiving synapse
                            Neuron[synConN].corticalStrength[l] = Neuron[synConN].corticalStrength[l] + updateStrength*gmax_cortical;
                            
                            // from synConN --> l (ONLY TRIPLET LTP)
                            // says if neuron l (o2) and synConN (r1) fired previously (o2_past) then do LTP based on timing between neuron foo spike (r1) and current tsp
                            updateStrength = Neuron[l].o2_past*Neuron[synConN].r1*A_trip_LTP_cort;
                            Neuron[l].corticalStrength[synConN] = Neuron[l].corticalStrength[synConN] + updateStrength*gmax_cortical;
                            // if too small, set to 0
                            if (Neuron[synConN].corticalStrength[l] < 0.0)
                            {Neuron[synConN].corticalStrength[l] = 0.0;}
                            // if too big, set to max
                            if (Neuron[l].corticalStrength[synConN] > gmax_cortical)
                            {Neuron[l].corticalStrength[synConN] = gmax_cortical;}
                        } // end if connected cell is exc
                        else // if connected cell is inhib, update strength according to that neuron's tracer
                        {
                            updateStrength_iSTDP = A_iSTDP*Neuron[synConN].iSTDP_trace; //convert tau from ms to sec
                            // update from synConN  (inhib) to l (exc)
                            Neuron[l].corticalStrength[synConN] = Neuron[l].corticalStrength[synConN] + updateStrength_iSTDP*gmax_cortical;
                            
                            if (Neuron[l].corticalStrength[synConN] > 2.0*gmax_cortical)
                            {  Neuron[l].corticalStrength[synConN] = 2.0*gmax_cortical;}
                        }
                    } // end loop over cortical connections
                }
            } // end if spiking neuron is exc.
            
            // if spiking neuron is INHIB
            else
            {
                Neuron[l].countInhibCorticalSpikes = Neuron[l].countInhibCorticalSpikes + 1.0; //update spike counter
                // loop over synaptic connections
                for (int newI=0; newI<numSynConn; newI++)
                {
                    // get neuron number of connected cell:
                    synConNI = Neuron[l].corticalConn[newI]; //neuron number of iith connection to neuron l
                    Neuron[synConNI].g_inhib1 = Neuron[synConNI].g_inhib1 + Neuron[synConNI].corticalStrength[l]*exp(-(t - Tsp)*oversigmaI); //update INHIB conductance
                    // STDP only if time is right
                    if (cortexLearnFlag)
                    {
                        // iSTDP if connected cell is exc
                        if (Neuron[synConNI].type == 0)
                        {
                            // now update plasticity from I --> E (l --> synConN)
                            updateStrength = A_iSTDP*(Neuron[synConNI].iSTDP_trace - 1.5*targetRate*(tau_inhib/1000.0)); //convert tau from ms to sec
                       //     cout << "update strength = " << updateStrength << endl;
                            // strength from l --> synConN
                            Neuron[synConNI].corticalStrength[l] = Neuron[synConNI].corticalStrength[l] + updateStrength*gmax_cortical;
                            // if too big, set to max
                            if (Neuron[synConNI].corticalStrength[l] > 2.0*gmax_cortical)
                            {Neuron[synConNI].corticalStrength[l] = 2.0*gmax_cortical;}
                            if (Neuron[synConNI].corticalStrength[l]<0.0)
                            {Neuron[synConNI].corticalStrength[l] = 0.0;}
                        } // end if connected cell is exc
                    }
                }
            }
            
            // how many GJ connections?
            numGJ = Neuron[l].GJConn.size();
            // if have GJ connection --
            // make sure to add the spikelet and update plasticity if THAT neuron went over threshold
            if (numGJ > 0)
            {
                GJcell = Neuron[l].GJConn[0]; //its GJ connected neuron... ONLY ONE IN THIS SIM
                if (Neuron[GJcell].v > vReset) // if this neuron was NOT just reset
                    {Neuron[GJcell].v = Neuron[GJcell].v + spikeletSize;}// add 1 mV spikelet
                    if (Neuron[GJcell].v > vT_local) // if this puts Neuron GJcell over threshold, reset and record spike time
                    {   // update tracer variables
                        Neuron[GJcell].o1 = Neuron[GJcell].o1 + 1.0;
                        Neuron[GJcell].o2 = Neuron[GJcell].o2 + 1.0;
                        Neuron[GJcell].r1 =  Neuron[GJcell].r1 + 1.0;
                        Neuron[GJcell].iSTDP_trace =  Neuron[GJcell].iSTDP_trace + 1.0;

                        // reset voltage
                        Neuron[GJcell].v = vReset; // set this neuron's voltage to reset.
                        Neuron[GJcell].tspN.push_back(t+dt); // put the end of the time step as the spike time
                        // update LGN --> cortex synapses
                        numLGNconnect = Neuron[GJcell].synapticConn.size();
                        for (int L = 0; L<numLGNconnect; L++)
                        {
                            // LGN index that is connected to neuron KK
                            indexL = Neuron[GJcell].synapticConn[L];
                            // length of synaptic spikes (should always be less than 5)
                            lengthTsp = Neuron[GJcell].tsp_input[indexL].size();
                            //  if LGN spiked previously
                            if (lengthTsp > 0)
                            {   //time difference from pre to post
                                    updateStrength = r1[indexL]*A_trip_LTP*Neuron[GJcell].o2_past;
                                    Neuron[GJcell].synapStrength[indexL] = Neuron[GJcell].synapStrength[indexL] + updateStrength*gmax;
                            }
                            // if weight goes above max, set to max:
                            if (Neuron[GJcell].synapStrength[indexL] > gmax)
                            { Neuron[GJcell].synapStrength[indexL] = gmax;}
                        }
                        // update synaptic connections for neuron GJcell
                        numSynConn_elect = Neuron[GJcell].corticalConn.size(); // get number of synaptic connections
                        // loop over thse post-synaptic cells:
                        for (int ii=0; ii<numSynConn_elect; ii++)
                        {
                            // get neuron number of connected cell:
                            foo_synConN = Neuron[GJcell].corticalConn[ii]; //neuron number of iith connection to neuron l
                            Neuron[foo_synConN].g_excite1 = Neuron[foo_synConN].g_excite1 + Neuron[foo_synConN].corticalStrength[GJcell]*exp(-(t - Tsp)*oversigmaE); //update EXC conductance
                            // STDP only if time is right
                            if (cortexLearnFlag)
                            {  //if the connected neuron is exc
                                if (Neuron[foo_synConN].type == 0)
                                {
                                    //  plasticity from GJcell --> foo_synConN (ONLY LTD)
                                    updateStrength = -Neuron[foo_synConN].o1*Neuron[foo_synConN].LTD_cort; // LTD should be for the neuron who is receiving synapse
                                    Neuron[foo_synConN].corticalStrength[GJcell] = Neuron[foo_synConN].corticalStrength[GJcell] + updateStrength*gmax_cortical;
                                  
                                    updateStrength2 = Neuron[GJcell].o2_past*Neuron[foo_synConN].r1*A_trip_LTP_cort;
                                    Neuron[GJcell].corticalStrength[foo_synConN] = Neuron[GJcell].corticalStrength[foo_synConN] + updateStrength2*gmax_cortical;
                                   
                                    // if too small, set to 0
                                    if (Neuron[foo_synConN].corticalStrength[GJcell] < 0.0)
                                    {Neuron[foo_synConN].corticalStrength[GJcell] = 0.0;}
                                    // if too big, set to max
                                    if (Neuron[GJcell].corticalStrength[foo_synConN] > gmax_cortical)
                                    {Neuron[GJcell].corticalStrength[foo_synConN] = gmax_cortical;}
                                } // end if connected cell is exc
                                else // if connected cell is inhib, update strength according to that neuron's tracer
                                {
                                    updateStrength_iSTDP = A_iSTDP*Neuron[foo_synConN].iSTDP_trace; //convert tau from ms to sec
                                    // update from synConN  (inhib) to l (exc)
                                    Neuron[GJcell].corticalStrength[foo_synConN] = Neuron[GJcell].corticalStrength[foo_synConN] + updateStrength_iSTDP*gmax_cortical;
                                    
                                    if (Neuron[GJcell].corticalStrength[foo_synConN] > 2.0*gmax_cortical)
                                    {  Neuron[GJcell].corticalStrength[foo_synConN] = 2.0*gmax_cortical;}
                                }
                            }
                        }
                    } // end if neuron foo spiked
            } // end if GJ coupled
            Neuron[l].Espike = false;
        } // end if spike
    } // end for loop over neurons
}

void poissonBackground(Network* Neuron, int N, double t, double dt, double oversigmaE, double backgroundNu)
{
    double t_SP = t;
    double r, u;
    int n;
    while (t_SP <= t+dt)
    {
        r = rand();
        u = r/RAND_MAX;
        while (u == 0)
        {
            r = rand();
            u = r/RAND_MAX;
        }
        t_SP = -log(u)/(N*backgroundNu) + t_SP; //Poisson spike time
        
        if (t_SP <= t+dt)
        {
            n  = rand()%N; //which neuron will recieve the spike
            Neuron[n].g_excite1 = Neuron[n].g_excite1 + Neuron[n].backgroundF*exp(-(t - t_SP)*oversigmaE);
            Neuron[n].countExternalSpikes = Neuron[n].countExternalSpikes + 1;
        }
    }
}

#endif /* learningNetworkModel_hpp */
