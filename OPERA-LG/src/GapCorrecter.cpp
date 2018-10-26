#include "GapCorrecter.h"

GapCorrecter::GapCorrecter(void){
	m_gapList = new vector<Gap*>;
}

GapCorrecter::GapCorrecter( int max, int min, int contigLength, int mean, int std, int readLength, vector<double> *pos ){
	m_gapList = new vector<Gap*>();
	m_contigLength = 0;
	m_contigLength = contigLength;
	m_maxGapSize = max;   // should be the max insert size - read length
	m_minGapSize = 0;
	m_minGapSize = min;
	m_mean = mean;
	m_std = std;

	m_possibility = pos;

	m_readLength = 0;
	m_readLength = readLength;
	m_useNormalDistribution = true;

	int sum = 0;
	for( int i = 0; i < (int) pos->size(); i++ ){
		sum += pos->at( i );
	}
	if( sum <= 100000 ){
		//cerr<<"use normal distribution\n";
		m_useNormalDistribution = true;
	}
	else
		m_useNormalDistribution = false;

	Init();
}


GapCorrecter::~GapCorrecter(void){
	for( int i = 0; i < (int) m_gapList->size(); i++ ){
		delete m_gapList->at( i );
	}
	m_gapList->clear();
	delete m_gapList;
}

void GapCorrecter::Init(){
	for( int i = 0; i <= m_maxGapSize; i++ ){
		Gap *newGap = new Gap();
		newGap->m_realGapSize = i; 
		m_gapList->push_back( newGap );
	}
}

// construct the relationship between real gap size and estimated gap size
void GapCorrecter::CalculateTable(){
	//cerr<<"mean: "<<m_mean<<endl;
	m_minEstimatedGapSize = 0;
	m_maxEstimatedGapSize = 0;

	//bool startIsZero = true;
	//bool endIsZero = false;

	// for gap = 0, 
	   // calculate summation of i*P(i) from g to g+c
	   // calculate summation of P(i) from g to g+c
	double numerator = 0;
	double denominator = 0;
	for( int i = m_readLength; i <= m_contigLength; i++ ){
		double possibility = CalculatePossibility( i );
		//double possibility = m_possibility->at( i );
		numerator += i * possibility;
		denominator += possibility;
	}

	   // calculate new mu
	double newMu;
	if( denominator == 0 )
		newMu = 0;
	else
		newMu = numerator / denominator;

	   // calculate estimated gap size
	   // g_head = g - mu_head + mu
	double newGapSize;
	if( m_contigLength >= m_minGapSize )
		newGapSize = 0 - newMu + m_mean;
	else
		newGapSize = 0;
	m_gapList->at( 0 )->m_estimatedGapSize = newGapSize; 
	m_gapList->at( 0 )->m_newMu = newMu;

	if( newGapSize != 0 )
		m_minEstimatedGapSize = newGapSize;

	if( newGapSize > m_maxEstimatedGapSize ) 
		m_maxEstimatedGapSize = newGapSize;

	// for gap g from 1 to max
	for( int i = 1; i < (int) m_gapList->size(); i++ ){
		double possibilityGM1 = CalculatePossibility( i - 1 + m_readLength); 
		double possibilityGPC = CalculatePossibility( i + m_contigLength );
		//double possibilityGM1 = m_possibility->at( i - 1 );
		//double possibilityGPC = m_possibility->at( i + m_contigLength );
		
		// update i*P(i): minus (g-1)P(g-1), plus (g+c)P(g+c)
		numerator -= (i - 1 + m_readLength) * possibilityGM1; 
		numerator += (i + m_contigLength) * possibilityGPC; 
		
		// update P(i): minus P(g-1), plus P(g+c)
		denominator -= possibilityGM1;
		denominator += possibilityGPC;

		// calculate new mu
		if( denominator > 0 )
			newMu = numerator / denominator;
		else
			newMu = 0;

		// calculate estimated gap size
		double newGapSize;
		if( denominator > 0 )
			newGapSize = i - newMu + m_mean;
		else
			newGapSize = 0;

		if( newGapSize > m_maxEstimatedGapSize ) 
			m_maxEstimatedGapSize = newGapSize;

		if( m_minEstimatedGapSize == 0 && newGapSize != 0 )
			m_minEstimatedGapSize = newGapSize;

		m_gapList->at( i )->m_estimatedGapSize = newGapSize; 
		m_gapList->at( i )->m_newMu = newMu;
	}
}

// calculate the possibility of a certain value
double GapCorrecter::CalculatePossibility( int value ){
	if( m_useNormalDistribution ){
		// use normal distribution
		double result = exp( -pow(value - m_mean, 2)/(2 * m_std * m_std) );
		return result * 10000;
	}
	else{
		double result;
		if( value >= (int) m_possibility->size() )
			result = 0;
		else
			result = m_possibility->at( value );
		
		return result;
	}
}

// print the constructed table
int GapCorrecter::PrintTable( string fileName ){
	ofstream tableWriter(  (Configure::OUTPUT_FOLDER + fileName).c_str() );
	if( tableWriter == NULL ){
		cout<<"ERROR: Cannot open "<<(Configure::OUTPUT_FOLDER + fileName)<<" file"<<endl;
		return -1;
	}
	
	string results;
	results = "read_gap_size\testimated_gap_size\tnew_mean\n";
	for( int i = 0; i < (int) m_gapList->size(); i++ ){
		results.append( itos( i ) + "\t" + itos( m_gapList->at( i )->m_estimatedGapSize ) + 
				"\t" + itos( m_gapList->at( i )->m_newMu ) + "\n" );
	}
	
	tableWriter.write( results.c_str(), results.length() );
	tableWriter.close();
	
	return 1;
}              


// get the corresponding real gap size of a given estimated gap size   
int GapCorrecter::GetRealGapSize( double estimatedGapSize )
{
	//cout<<"Estimated gap size is: "<<estimatedGapSize<<endl;

	//cerr<<"Get real gaps size\n";
	// if the value is too small
	if( estimatedGapSize < m_minEstimatedGapSize )
		return estimatedGapSize;

	// if the value is too big
	if( estimatedGapSize > m_maxEstimatedGapSize ){
		// bigger than the tail value, just use the tail value
		//int number = 0;
		//double sum = 0;
		
		for( int i = m_gapList->size() - 1; i >= 0; i-- ){
			if( i == 0 ){
				continue;
			}

			if( m_gapList->at( i )->m_estimatedGapSize > 0 ){ //== m_maxEstimatedGapSize ){
				return i;
				//number++;
				//sum += i;
			}
			//else
			//	break;
		}
		//double realValue = sum / number;
		//return realValue;
		return estimatedGapSize;
	}

	// search for the right position
	int max = m_gapList->size() - 1;
	int min = 0;

	int middle = -1;

	while( min != max ){
		middle = (min + max)/2;
		if( m_gapList->at( middle )->m_estimatedGapSize == estimatedGapSize ){
			// find the exact value, exit
			break;
		}

		if( middle == min || middle == max )
			break;
		
		if( m_gapList->at( middle )->m_estimatedGapSize < estimatedGapSize ){
			min = middle;
		}
		else{
			max = middle;
		}
	}

	middle = (min+max)/2;

	// collect all the possible gap sizes
	int number = 1;
	double sum = middle;
    
	// check the one before middle
	for( int i = middle - 1; i >= 0; i-- ){
		if( m_gapList->at( i )->m_estimatedGapSize == m_gapList->at( middle )->m_estimatedGapSize ){
			number++;
			sum += i;
		}
		else
			break;
	}

	// check the one after middle
	for( int i = middle + 1; i < (int) m_gapList->size(); i++ ){
		if( m_gapList->at( i )->m_estimatedGapSize == m_gapList->at( middle )->m_estimatedGapSize ){
			number++;
			sum += i;
		}
		else
			break;
	}
	
	//cerr<<"find the value\n";
	double realValue = sum / number;
	return realValue;
}     


