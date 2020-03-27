//
// tripletRule_fullSim_iSTDP_diffates_radConn
//
// triplet rule with homeostatic drive
// includes inhibitory neurons and iSTDP from I --> E
// Radius of connectivity for the exc. connectivity
//
//
//  Created by Jen Crodelle on 6/10/2019
//
//
#include "trip_iSTDP_funcs.hpp"
#include <iostream>

int main(int argc, char * argv[]){
    // input params
    double seed = atof(argv[1]);    //random seed
    double Tfin = atof(argv[2]);  // final time (ms)
    int N = atof(argv[3]); // number of neurons
    double backgroundF =atof(argv[4]); //background strength
    double extRate =atof(argv[5]); //external rate
    double input_A_tripLTP =atof(argv[6]); // triplet LTP
    double A_trip_LTP_cort =atof(argv[7]); // triplet LTP
    double A_iSTDP =atof(argv[8]); // iSTDP weight
    double gmax =atof(argv[9]);  // max weight strength
    double probGJs = atof(argv[10]); // prob of connecting with GJ
   // double probGJs = atof(argv[11]);//
    double targetRate = atof(argv[11]); // target rate for changing LTD
    double gmax_cortical = atof(argv[12]); // max cortical weight
    double timeForSynapses =atof(argv[13]); //time to start recurrent synapses
    
    if (argc != 14)
    { cout << " there are " << argc << " input arguments, but there should be 14!" << endl;
        return 0; }
// setting up simulation params
    srand(seed);
    int T = (int)20*Tfin;
    double dt = Tfin/T; //step size
    cout << "final time  = " << Tfin << endl;
    cout << "dt = " << dt << " T  = " << T << endl;
    
    int radiusConn[2][2] = {{4,10},{10,10}};
    double probCortConnect = 1.0; //prob of cortical cells connecting
    //double strengthConnect = 0.0; // strength of initial connection
//    int radiusConn[2][2] = {{40,126},{126,126}};
    // E-E,
    double strengthConnect[2][2] = {{0.0*gmax_cortical,0.1*gmax_cortical},{0.3*gmax_cortical,0.1*gmax_cortical}};
    Network* Neuron = new Network[N];

    // create neurons
    create_neurons(Neuron, N, T, gmax, N_input, backgroundF);
    // create neurons with LGN seeded
    //create_neurons_seedLGN(Neuron, N, T, gmax, N_input, backgroundF);
    cout << "there are " << N << " neurons with " << N_input << " synapses" << endl;
   
    // just set up the syapses as all-to-all with probability probCortConnect
    setUpZeroSynapses(Neuron, N,  probCortConnect);
    cout << " all-to-all synapses " << endl;
    parameterFile << "all-to-all synapses" << endl;
    // set up synaptic structure with radius
   // setUpZeroSynapses_radius(Neuron,N,radiusConn);
   //  cout << " radius synapses " << endl;
  //  parameterFile << "radius synapses" << endl;
    //double strengthConnect[2][2] = {{0.0,0.1*gmax_cortical},{0.25*gmax_cortical,0.1*gmax_cortical}};
    //double probGJs = 0.0; //0.02; // prob of connecting to sister cell
    // connect with gap junctions
    int numGJs = 150;
    //connectGJs(Neuron, N,  numGJs);
    //parameterFile << "pairs of GJs" << endl;
    
    connectGJs_sisters(Neuron, N, probGJs);
    parameterFile << "sister GJs" << endl;
    // set initial interval length:
    double  r = rand();
    double  u = r/RAND_MAX;
    bool errorflag = false;
    ::t_interval = -log(u)*(20.0); //initial interval..
    
    double countForIntervals = 0; //count number of times we update the interval length
    bool cortexLearnFlag;  // flag for when to turn cortical learning on
    double updateStrength, numGJconn, sumGJ, gc, spikeletSize,vT_local,A_trip_LTP;
    int indexL, numLGNconnect, lengthTsp, foo;
    
    int onceFlag = 1;
    //i increments with timestep
    for (int i = 0; i < T; i++)
    {
        if (i%200000 == 0) //only record sometimes the weights
        {
            cout << "time: " << t/1000.0 << "sec"<< endl;
            // record the weights for the LGN synapses and cortical neurons
            for (int KK = 0; KK < N; KK++)
            {
                for (int L = 0; L < N_input; L++)
                {
                    if (L==N_input-1)
                    {
                        // if(i==0) {initialWeights << Neuron[KK].synapStrength[L]<< endl;}
                        W << Neuron[KK].synapStrength[L] << endl;}
                    else
                    {//if(i==0){initialWeights << Neuron[KK].synapStrength[L]<< ",";}
                        W << Neuron[KK].synapStrength[L] << ",";}
                }
                for (int J = 0; J < N; J++)
                {
                    if (J==N-1)
                    { W_cortical << Neuron[KK].corticalStrength[J]  << endl;}
                    else { W_cortical << Neuron[KK].corticalStrength[J] << ",";}
                }
            }
        }
        // only allow neurons to learn after 3 seconds of background:
        if(t < 3000.0)
        {
            A_trip_LTP = 0.0;
            cortexLearnFlag = false;
           // gc = 0.004;
             gc = 0.06;
            spikeletSize = 1.0;
        }
        // only allow feedforward synapses and GJs
        else if (t < timeForSynapses)
        {
            A_trip_LTP = input_A_tripLTP;
            cortexLearnFlag = false;
          //  gc = 0.004;
            gc = 0.06;
            spikeletSize = 1.0;
//            if (onceFlag == 1)
//            {
//                connectSynapses(Neuron, N, strengthConnect,gmax_cortical);
//                onceFlag = 0;
//            }
        }
        // now feedforward AND recurrent NO GJs
        else
        {
            A_trip_LTP = input_A_tripLTP;
            cortexLearnFlag = true;
            gc = 0.0;
            spikeletSize = 0.0;
            // once we allow recurrent connections set synaptic connections
//            if (onceFlag == 1)
//            {
//                connectSynapses(Neuron, N, strengthConnect,gmax_cortical);
//                cout << "connected synapses at time t = " << t<< endl;
//                onceFlag = 0;
//            }
        }
        
        for(int L = 0; L<N_input; L++)
        {
            r1[L] = r1[L]*exp(-dt/tau_LTP);
        }
        
        for (int KK = 0; KK<N; KK++)
        {
            // update conductance:
            Neuron[KK].g_excite1 = Neuron[KK].g_excite0*exp(-dt*oversigmaE);
            Neuron[KK].g_inhib1 = Neuron[KK].g_inhib0*exp(-dt*oversigmaI);
            Neuron[KK].avgVolt = Neuron[KK].avgVolt*exp(-dt*overTau_v);
            // update tracer:
            Neuron[KK].o1 = Neuron[KK].o1*exp(-dt/tau_LTD);
            Neuron[KK].o2 = Neuron[KK].o2*exp(-dt/trip_tau_LTD);
            Neuron[KK].r1 = Neuron[KK].r1*exp(-dt/tau_LTP);
            Neuron[KK].iSTDP_trace = Neuron[KK].iSTDP_trace*exp(-dt/tau_inhib);
            
            // calculate GJ stuff:
            numGJconn = Neuron[KK].GJConn.size();
            if (numGJconn > 0)
            { // loop over GJ connections and calculate sum
                vT_local = -45.0; //threshold voltage
                sumGJ = 0.0;
                for (int J = 0; J<numGJconn; J++)
                {
                    foo = Neuron[KK].GJConn[J];
                    sumGJ = sumGJ + Neuron[foo].v;
                }
            }
            else
            {
                gc = 0.0;
                vT_local = -45.0; //threshold voltage
                sumGJ = 0.0;
                numGJconn = 0.0;
            }
            //calculate variables: a0,a1,b0,b1,k1,k2 for RK2 Method
            oldv = Neuron[KK].v;
            gE0 = Neuron[KK].g_excite0;
            gI0 = Neuron[KK].g_inhib0;
            
            a0 = (gL+gE0+gI0+gc*numGJconn)*overC;
            b0 = (gL*vR+gE0*vE+gI0*vI+gc*sumGJ)*overC;
            
            gE1 = Neuron[KK].g_excite1;
            gI1 = Neuron[KK].g_inhib1;
            
            a1 = (gL+gE1+gI1+gc*numGJconn)*overC;
            b1 = (gL*vR+gE1*vE+gI1*vI+gc*sumGJ)*overC;
            
            k1 = -a0*oldv+b0;
            k2 = -a1*(oldv+dt*k1)+b1;
            
            Neuron[KK].v = oldv+(dt/2)*(k1+k2); //RK2 step
            
            // check if voltage went below rest...
            if (Neuron[KK].v < (vR-5))
            {   cout << " We have a problem, v = " << Neuron[KK].v << endl;
                errorflag = true;
                break;
//                cout << " Neuron " << KK << " is type " << Neuron[KK].type << endl;
//                cout << " inhib conduct = " << Neuron[KK].g_inhib1 << endl;
//                cout << " old v = " << oldv << ", avg volt = " << Neuron[KK].avgVolt << endl;
            }
            // if reached threshold, neuron fired, update conductance and synaptic weights:
            if (Neuron[KK].v > vT_local)
            {
                // assume linear b/t vn and vn+1, solve for tspike:
                tspike = t+dt*((vT_local-oldv)/(Neuron[KK].v-oldv));
                // linearly interpolate to find new vn+1 value after spike - set back to -60 not -70:
                vtilda = (vReset - 0.5*(tspike-t)*(b0+b1-dt*a1*b0))/(1.0+0.5*(tspike-t)*(-a0-a1+dt*a0*a1));
                // re-evaute RK2 step the new voltage:
                k1tilda = -a0*vtilda+b0;
                k2tilda = -a1*(vtilda+k1tilda*dt)+b1;
                
                //replace with new spiked voltage value
                Neuron[KK].v = vtilda+(dt/2)*(k1tilda+k2tilda);
                
                // if this voltage is still over threshold, there's an issue..
                if (Neuron[KK].v > vT_local)
                {   cout << "Two spikes occured in one time step. Vtilda was above threshold." << endl;
                    errorflag = true;
                    break;
                }
                
                //collect spike time for this neuron
                Neuron[KK].tspN.push_back(tspike);
                Neuron[KK].avgVolt = Neuron[KK].avgVolt + 1000.*overTau_v;
                Neuron[KK].Espike = true; //mark that THIS neuron spiked
                
                numLGNconnect = Neuron[KK].synapticConn.size();
                for (int L = 0; L<numLGNconnect; L++)
                {
                    // LGN index that is connected to neuron KK
                    indexL = Neuron[KK].synapticConn[L];
                    // length of synaptic spikes (should always be less than 5)
                    lengthTsp = Neuron[KK].tsp_input[indexL].size();
                    // if LGN spiked previously
                    if (lengthTsp > 0)
                    {   //time difference from pre to post
                        updateStrength = r1[indexL]*A_trip_LTP*Neuron[KK].o2_past;
                        Neuron[KK].synapStrength[indexL] = Neuron[KK].synapStrength[indexL] + updateStrength*gmax;
                        //  cout << "LTP: strength = " << updateStrength<<endl;
                    }
                    // if weight goes above max, set to max:
                    if (Neuron[KK].synapStrength[indexL] > gmax)
                    { Neuron[KK].synapStrength[indexL] = gmax;}
                } // end for loop over synapses
            } // end if over threshold
        } // end loop over cortical cells
        
        // model incoming spikes from LGN and update synapses onto cortical cells
        if(t < ::t_interval) //still in the same interval, keep firing rate constant.
        { inputFromLGN(Neuron, N, N_input, t, dt, oversigmaE, tau_LTD, gmax, 0, extRate);}
        else //calculate new interval and firing rates
        {
            countForIntervals = countForIntervals+1;
            inputFromLGN(Neuron,N, N_input, t, dt, oversigmaE, tau_LTD, gmax, 1, extRate);
        }
        // BACKGROUND DRIVE
        poissonBackground(Neuron, N, t, dt, oversigmaE, backgroundNu);
        spikeUpdate(Neuron, N, t, dt, gmax, A_trip_LTP_cort, A_trip_LTP, cortexLearnFlag, A_iSTDP, tau_inhib, gmax_cortical, spikeletSize,targetRate);
        
        //shift all current conductance values to be old g calculate new g at the next timestep
        for (int c = 0; c < N; c++)
        {
            Neuron[c].g_excite0 = Neuron[c].g_excite1; //shift conductance values
            Neuron[c].g_inhib0 = Neuron[c].g_inhib1; //shift conductance values
            Neuron[c].o2_past = Neuron[c].o2; // update past o2 trace
            // update the LTD for each postsynaptic neuron
            Neuron[c].LTD = Neuron[c].avgVolt*Neuron[c].avgVolt*(A_trip_LTP*tau_LTP*trip_tau_LTD)/(targetRate*1000.0*tau_LTD);
            Neuron[c].LTD_cort =  Neuron[c].avgVolt*Neuron[c].avgVolt*(A_trip_LTP_cort*tau_LTP*trip_tau_LTD)/(targetRate*1000.0*tau_LTD);
        }
        if (errorflag == true)
        {break;}
        //Collect Data from this timestep
        t = t+dt;
    }
    // end time loop.
    ////
    ////
    ////
    ////
    
    // read in spike times for each neuron:
    double avgInputSpikes = 0.0;
    double avgExternalSpikes = 0.0;
    double avgExcCorticalSpikes = 0.0;
    double avgInhibCorticalSpikes = 0.0;
    for (int k = 0; k<N; k++)
    {
        // count cortical spikes and LGN spikes
        avgInputSpikes = avgInputSpikes + Neuron[k].countInputSpikes;
        avgExternalSpikes = avgExternalSpikes + Neuron[k].countExternalSpikes;
        double lengthSp = Neuron[k].tspN.size();
        if (Neuron[k].type == 0)
        {avgExcCorticalSpikes = avgExcCorticalSpikes + lengthSp;}
        else{avgInhibCorticalSpikes = avgInhibCorticalSpikes + lengthSp; }
        // enter neurons spike times into file:
        spTimes << Neuron[k].type << "," <<  Neuron[k].sisterID << "," << lengthSp << ",";
        if (lengthSp == 0)
        { spTimes << endl;
        }
        for (int j = 0; j<lengthSp; j++)
        {
            if (j == lengthSp - 1)
            { spTimes << Neuron[k].tspN[j] << endl;
            }
            else
            { spTimes << Neuron[k].tspN[j] << ","; }
        }
        // record the numer of electric connections
        int numElect = Neuron[k].GJConn.size();
        electConn_file << numElect << ",";
        if (numElect > 0)
        {
        for (int c=0; c<numElect; c++)
        {
            if (c == numElect-1)
            { electConn_file << Neuron[k].GJConn[c] << endl;}
            else
            { electConn_file << Neuron[k].GJConn[c] << ",";}
        }
        }
        else
        {electConn_file << endl;}
        
        // for the synaptic strengths from LGN to cortex
        for (int L = 0; L<N_input; L++)
        {
            if(L == N_input-1)
            {   W << Neuron[k].synapStrength[L] << endl;
                finalWeights << Neuron[k].synapStrength[L] << endl;}
            else
            {    W << Neuron[k].synapStrength[L] << ",";
                finalWeights << Neuron[k].synapStrength[L] << ",";}
        }
    }
    
    avgInputSpikes = (1000.0*avgInputSpikes)/((double)N_input*probSynapseConn*N*Tfin);
    avgExcCorticalSpikes =(1000.0*avgExcCorticalSpikes)/((double)(N-numInhib)*Tfin);
    avgInhibCorticalSpikes =(1000.0*avgInhibCorticalSpikes)/((double)numInhib*Tfin);
    cout << " Avg incoming FR from LGN cells = " << avgInputSpikes << endl;
    cout << " Avg exc cortical FR  = " << avgExcCorticalSpikes << endl;
    cout << " Avg inhib cortical FR  = " << avgInhibCorticalSpikes << endl;

    // parameters in a specific order to be read into continuation file
    parameterFile << "N = " << N << endl;
    parameterFile<< "Final time = " << Tfin << endl;
    parameterFile<< "dt = " << dt << endl;
    parameterFile << "background F = " << backgroundF << endl;
    parameterFile << " background Nu = " << backgroundNu << endl;
    parameterFile << " external rate R0 = " << extRate << endl;
    parameterFile << " prob of GJ = " << probGJs << endl;
    parameterFile << "gc = 0.06" << endl;
    parameterFile << " rand seed = " << seed << endl;
    
    parameterFile << "radius E-E" <<radiusConn[0][0] << endl;
    parameterFile << "radius I-E" <<radiusConn[0][1] << endl;
    parameterFile << "radius E-I" <<radiusConn[1][0] << endl;
    parameterFile << "radius I-I" <<radiusConn[1][1] << endl;
    parameterFile<< "gmax_cortical = " << gmax_cortical << endl;
    parameterFile << "prob cortical connection  " << probCortConnect << endl;

    parameterFile<<" A3_LTP = " << A_trip_LTP << endl;
    parameterFile << "A3 LTP cort = " << A_trip_LTP_cort << endl;
    parameterFile<< "tau_LTD = " << tau_LTD << endl;
    parameterFile << "tau_LTP = " << tau_LTP << endl;
    parameterFile << "trip_tau_LTD " << trip_tau_LTD << endl;
    parameterFile<< "gmax = " << gmax << endl;
    //parameterFile << "number GJs = " << numGJs << endl;
    parameterFile << "target rate = " << targetRate << endl;
    
    parameterFile << "gL = " << gL << endl;
    parameterFile << "Vrest = "  << vR << endl;
    parameterFile<< "Vthreshold = " << vT << endl;
    parameterFile << "Sigma_E =  " << sigma_E << endl;
    parameterFile << "Sigma_E=I =  " << sigma_I << endl;
    
    parameterFile << "number of synapses " << N_input << endl;
    parameterFile<< "avg input from LGN neuron = " << avgInputSpikes << endl;
    parameterFile << " Avg exc cortical FR  = " << avgExcCorticalSpikes << endl;
    parameterFile << " Avg inhib cortical FR  = " << avgInhibCorticalSpikes << endl;

    parameterFile << "file run: tripletRule_fullSim_iSTDP_diffRates_radConn.cpp" << endl;


    W.close();
    finalWeights.close();
    spTimes.close();
    parameterFile.close();
    electConn_file.close();
    delete[] Neuron;
    
    cout << clock()/CLOCKS_PER_SEC << " seconds " << endl;
}

