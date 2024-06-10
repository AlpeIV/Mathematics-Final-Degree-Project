#include <stdio.h>
#include <ctime>
#include <stdlib.h>
#include <cmath>
#include <float.h>
#include <random>
#include <iostream>
#include <fstream>
#include <chrono>
#include <unistd.h>
#include <chrono>
#include <limits>
#include <algorithm>
#include <bits/stdc++.h> 





std::string RandomString(const int len) {
    static const char alphanumerics[62] = {
            0x30, 0x31, 0x32, 0x33, 0x34, 0x35, 0x36, 0x37, 0x38, 0x39, 0x41, 0x42,
            0x43, 0x44, 0x45, 0x46, 0x47, 0x48, 0x49, 0x4A, 0x4B, 0x4C, 0x4D, 0x4E,
            0x4F, 0x50, 0x51, 0x52, 0x53, 0x54, 0x55, 0x56, 0x57, 0x58, 0x59, 0x5A,
            0x61, 0x62, 0x63, 0x64, 0x65, 0x66, 0x67, 0x68, 0x69, 0x6A, 0x6B, 0x6C,
            0x6D, 0x6E, 0x6F, 0x70, 0x71, 0x72, 0x73, 0x74, 0x75, 0x76, 0x77, 0x78,
            0x79, 0x7A
    };
    srand(time(0));

    std::string tmp;
    tmp.reserve(len);
    for (int i = 0; i < len; ++i) {
        tmp += alphanumerics[rand() % (sizeof(alphanumerics) - 1)];
    }
    return tmp;
}





long double PowerDistribution(long double tailpower, long double SightLimit){
	

	std::mt19937_64 generator(rand());
	std::uniform_real_distribution<long double> distribution(0.0L,1.0L); //Standard Uniform
	
	long double UniformRandNumb = distribution(generator);
	long double x = (SightLimit * powl(1.0L - UniformRandNumb, -1.0L / (tailpower - 1.0L)));
	
	
	//std::cout << "Dist: " << x << std::endl;
	return x;
}


int LeftOrRight(){
	int Cointoss = rand() % 2;
	if(Cointoss){
		return Cointoss;
	}else{
		return -1;	
	}
	return 1;
}


long double UniFlight(long double tailpower, long double SightLimit){
	
	long double NewDispl = LeftOrRight() * PowerDistribution(tailpower, SightLimit);
	return NewDispl;
}



void ApplyLinearFlight(long double &x_axis, long double dist){
	x_axis = x_axis + dist;
	return;
}




long long NaturalToInteger(int Value){
	long long Resolt = Value % 2;
	if(Resolt){
		return -((Value + 1) / 2);
	}else{
		return Value/2;	
	}
	return 1;
}


long double LinearDistance(long double pos1, long double pos2){
	long double Resolt = pos1 - pos2;
	if(Resolt < 0){
		return -1.0L * Resolt * 1.0L;
	}
	return Resolt;
}


long double GlobalDistance = 0.0L, GlobalTries = 0.0L;




void Pseudomain(void){
	srand (time(NULL));
	
	std::vector<long double> NodePositions; //
	std::vector<long double> NodeDistances; //We assume a node at 0




	NodePositions.reserve(3);
    NodePositions.push_back(-100.0L);
    NodePositions.push_back(100.0L);

  	
  	
  	std::cout << "Sorting Positions..." << std::endl;
  	std::sort(NodePositions.begin(), NodePositions.end(),  std::less<long double>());
  	
	
	for (size_t i = 0; i < NodePositions.size(); ++i) {
        std::cout << "Pos: " << NodePositions[i] << " " << std::endl;
        
    }
    
    
      	 std::cout  << std::endl;
    
    
    long double SightRange = 0.25L;   //r_v
    long double tailPower = 2.0L;     //Manual Input
    
    std::ofstream fout;
	std::string FileName;
	if(tailPower > 3.0){
		FileName = "NonDestructiveBrownian1DWalk_" + RandomString(16);
	}else{
		FileName = "NonDestructiveLevyMotion_" + RandomString(16);
	}
	
	fout.open(FileName.c_str());
	
    
    
	
	std::cout << "Sight Range chosen (r_v): " << SightRange << " Distance Units" << std::endl;
	std::cout << "Power Law Tail: " << tailPower << "" << std::endl;
	
	int MaxIter = 10000;
  	std::cout << "Computing " << MaxIter << " Forager Flight Simulations"  << std::endl;
  	
  	int CurrentIter = 0;
  	long long TotalFlights = 0;
  	long double traversedDistance = 0.0L;
  
  
  	 std::cout  << std::endl;

  	
  	do{
  		
  		long double foragerPosition = 0.0L;
  		//std::cout << "Forager Starting Position: " << foragerPosition << " " << std::endl;
  		
  		long long LocalFlights = 0;
  		long double LocaltraversedDistance = 0.0L;
  		
  		int i = 0;
  		
		size_t Index = i;
  		
  		//Force a first flight from the collected node at origin before restoring:
  		
  			  	LocalFlights++;
  				long double NewDispl2 = UniFlight(tailPower , SightRange);
  				long double FuturePosition2 = foragerPosition + NewDispl2;
  				
  				//std::cout << foragerPosition << " forager pos" << std::endl;
  				//std::cout << FuturePosition2 << " next forager pos if no nodes" << std::endl;
  			
  				if((FuturePosition2 > (NodePositions[Index + 1] - SightRange))){
  					//we went beyond a node, so forager will stop at the node itself
  					LocaltraversedDistance = LocaltraversedDistance + (NodePositions[Index + 1] - foragerPosition);
  					goto Donee;
				}else if((FuturePosition2 < (NodePositions[Index] + SightRange))){
					LocaltraversedDistance = LocaltraversedDistance + (foragerPosition - NodePositions[Index]);
					goto Donee;
  					break;
				}else{
					LocaltraversedDistance = LocaltraversedDistance + fabs(NewDispl2);
  					foragerPosition = FuturePosition2;
				}
  				
		
  		
  		//Restore origin
  		NodePositions.push_back(0.0L);
  		std::sort(NodePositions.begin(), NodePositions.end(),  std::less<long double>());
  		
  		
  		//Select closer nodes at start position
  		i = -1;
  		do{
  			i++;
		}while(!((foragerPosition > NodePositions[i]) && (foragerPosition < NodePositions[i+1])));
		Index = i;
  		
  		
  		//std::cout << "Closer Nodes: " << NodePositions[i] << " & " << NodePositions[i + 1] << std::endl;
  	
  	
  	    	
  			while((foragerPosition > NodePositions[Index] + SightRange) && (foragerPosition < NodePositions[Index + 1] + SightRange)){ //We inbetween
  				LocalFlights++;
  				long double NewDispl = UniFlight(tailPower , SightRange);
  				long double FuturePosition = foragerPosition + NewDispl;
  			
  				//std::cout << foragerPosition << " forager pos" << std::endl;
  				//std::cout << FuturePosition << " next forager pos if no nodes" << std::endl;
  			
  				if((FuturePosition > (NodePositions[Index + 1] - SightRange))){
  					//we went beyond a node, so forager will stop at the node itself
  					LocaltraversedDistance = LocaltraversedDistance + (NodePositions[Index + 1] - foragerPosition);
  					break;
				}else if((FuturePosition < (NodePositions[Index] + SightRange))){
					LocaltraversedDistance = LocaltraversedDistance + (foragerPosition - NodePositions[Index]);
  					break;
				}else{
					LocaltraversedDistance = LocaltraversedDistance + fabs(NewDispl);
  					foragerPosition = FuturePosition;
				}
			};
		
		Donee:
		
			//std::cout << "Took: " << LocalFlights << " flights" << std::endl;
  			TotalFlights = TotalFlights + LocalFlights;
  			//std::cout << "Traversed: " << LocaltraversedDistance << " Distance" << std::endl;
  			traversedDistance = traversedDistance + LocaltraversedDistance;
  			//std::cout  << std::endl;	
			CurrentIter++;	
			
			if((CurrentIter%100) == 50){  //Scoped, save results every some successes
						fout << std::endl;
						fout << std::endl;
						fout << "Success Backup: " << std::endl;
						fout << TotalFlights << " flights" << std::endl;
						fout << CurrentIter << " successes" << std::endl;
						fout << "Traversed " << traversedDistance << " Distance Units" << std::endl;
  	
  						double AverageFlightCount = (1.0 * TotalFlights * 1.0)/(1.0 * CurrentIter * 1.0);
  						long double AverageSuccessLength = traversedDistance/(1.0L * CurrentIter * 1.0L);
  						long double AverageFlightLength = traversedDistance/(1.0L * TotalFlights * 1.0L);
  	
  						fout << "Current Averages: " << std::endl;
						fout << AverageFlightCount << " Flights per success" << std::endl;
						fout << AverageSuccessLength << " Distance Units per success" << std::endl;
  						fout << AverageFlightLength << " Distance Units per flight" << std::endl;
						fout << std::endl;
						fout << std::endl;
			}
		
		
  		
	}while(CurrentIter < MaxIter);
	
	
	
		std::cout << "Sight Range chosen (r_v): " << SightRange << " Distance Units" << std::endl;
	std::cout << "Power Law Tail: " << tailPower << "" << std::endl;
	
	std::cout  << std::endl;
	std::cout << "For: " << CurrentIter << " tries" << std::endl;
	std::cout << "it took: " << TotalFlights << " flights" << std::endl;
	std::cout << "For a total of " << traversedDistance << " Distance Units Traversed" << std::endl;
	 std::cout  << std::endl;	
  	
  	double AverageFlightCount = (1.0 * TotalFlights * 1.0)/(1.0 * CurrentIter * 1.0);
  	long double AverageSuccessLength = traversedDistance/(1.0L * CurrentIter * 1.0L);
  	long double AverageFlightLength = traversedDistance/(1.0L * TotalFlights * 1.0L);
  	
  	std::cout << "For an average of: " << std::endl;
	std::cout << AverageFlightCount << " flights" << std::endl;
	std::cout << AverageSuccessLength << " Distance Units per success" << std::endl;
  	std::cout << AverageFlightLength << " Distance Units per flight" << std::endl;
  	
  	
  	std::cout << std::endl;
   	std::cout << std::endl;
	std::cout << std::endl;
	
	
							fout << std::endl;
						fout << std::endl;
												fout << std::endl;
						fout << std::endl;
	fout << std::endl;
						fout << std::endl;
						fout << "\\mu = " << tailPower << std::endl;
						fout << TotalFlights << " flights" << std::endl;
						fout << CurrentIter << " successes" << std::endl;
						fout << "Traversed " << traversedDistance << " Distance Units" << std::endl;
  	
  						AverageFlightCount = (1.0 * TotalFlights * 1.0)/(1.0 * CurrentIter * 1.0);
  						AverageSuccessLength = traversedDistance/(1.0L * CurrentIter * 1.0L);
  						AverageFlightLength = traversedDistance/(1.0L * TotalFlights * 1.0L);
  	
  						fout << "Current Averages: " << std::endl;
						fout << AverageFlightCount << " Flights per success" << std::endl;
						fout << AverageSuccessLength << " Distance Units per success" << std::endl;
  						fout << AverageFlightLength << " Distance Units per flight" << std::endl;
						fout << std::endl;
						fout << std::endl;
	
	
			 
			  	
  	return ;
  	
  	
  	
  	
  	
  	
  	
  	//Miscellaneous computations
  	
  	
  	std::cout << "Computing Distance between Nodes..." << std::endl;
  	for (size_t i = 0; i < NodePositions.size() - 1; ++i) {
        long double difference = NodePositions[i + 1] - NodePositions[i];
        NodeDistances.push_back(difference);
    }
    
    /*
    for (size_t i = 0; i < NodeDistances.size(); ++i) {
        std::cout << "Distances: " << NodeDistances[i] << " " << std::endl;
        
    }
    */

    //Average distance in a single region is (Nodes + 1)/Chosen area; if 0 nodes
  	long double AverageDistance = 0.0L;
	std::cout << "Computing Average Distance..." << std::endl;
  	for (size_t i = 0; i < NodeDistances.size(); ++i) {
        AverageDistance = AverageDistance + NodeDistances[i];
    }
    
  	AverageDistance = AverageDistance/(1.0L * NodeDistances.size() * 1.0L);
  	std::cout << "Computed Average Distance: " << AverageDistance << " Units" << std::endl;
	
	
    	
    	

	
	
	
	return ;
	
	long double x_axis = 0.0L;

	
	long double PowerTail = 2.0L;
	long double MaxSight = 1.0L;
	
	for(unsigned int i = 0; i < 30; i++){
		long double LineFlight = UniFlight(PowerTail, MaxSight);
		ApplyLinearFlight(x_axis, LineFlight);
		std::cout << "NewPos: " << x_axis << std::endl;
	}
	
	
	
	
	
	
	
	
	return ;
}





int main(void){

		Pseudomain();

	
	return 0;
}


