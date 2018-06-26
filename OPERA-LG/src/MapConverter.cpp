#include "MapConverter.h"

MapConverter::MapConverter(void)
{
	m_singlePetsMap = new map<pair<int, int>, multiset<SinglePet*, less_distance>*>;
	m_pairMap = new list<string>;
	m_totalTime = 0;

	m_libString = "";

	m_headHashMap = new hash_map<const char*, list<int>*, hash<const char*>, eqName>;
	m_tailHashMap = new hash_map<const char*, list<int>*, hash<const char*>, eqName>;
}

MapConverter::~MapConverter(void)
{
	delete m_singlePetsMap;
	delete m_pairMap;
	m_ClustersMap->clear();
	delete m_ClustersMap;

	for( int i = 0; i < (int)m_gapTables.size(); i++ ){
		delete m_gapTables.at( i );
	}
}


MapConverter::MapConverter( Graph *graph ){
	m_graph = graph;
	m_singlePetsMap = new map<pair<int, int>, multiset<SinglePet*, less_distance>*>;
	m_pairMap = new list<string>;
	m_totalTime = 0;
	m_pairTime = 0;
	
	m_ClustersMap = new map<pair<string, string>, multiset<PET*, less_std>*>;
}


// analyze mapping file or opera edge file
int MapConverter::AnalyzeMultiLib( string fileName, list<PetLibrary*> *libs  ){
	/*switch( Configure::MAP_TYPE ){
	case BOWTIE:
	case SAM:
		return AnalyzeBowtieMultiLib( libs );
	case OPERA:
		return AnalyzeOperaMultiLib( libs );
	default:
		return -1;
		}*/
	//Export the mapping based edges if any
	if( AnalyzeBowtieMultiLib( libs ) ); // To silence the return value
	//Add any OPERA edge file to libs
	AnalyzeOperaMultiLib( libs );
	return 1;
}


bool LibSort( PetLibrary *&lib1, PetLibrary *&lib2 ){
	return lib1->GetMean() < lib2->GetMean();
}

// analyze bowtie file
int MapConverter::AnalyzeBowtieMultiLib( list<PetLibrary*> *libs ){
	// get all mapping files and handle them one by one
	int i = 0;
	for( list<PetLibrary*>::iterator libIter = libs->begin(); libIter != libs->end(); libIter++ )
	{
		i++;
		cerr<<"Current library: "<<i<<" out of "<<Configure::MULTI_LIB_INFO->size()<<endl;

		if( ConvertBowtieFileMultiLib( *libIter ) == -1 )
			return -1;
	}
	
	// sort according to mean
	libs->sort( LibSort );
	return 1;
}

// read and analyze contig graph
// Useless function? ##?
int MapConverter::ConvertSoapContigGraph( string fileName, list<PetLibrary*> *libs )
{
	// usedContig, heads, tails
	// save the mappings of contig to scaffold
	//map<int, string> *heads;
	//map<int, string> *tails;

	// read contig file and save all scaffold and C id
	//list<Contig*> *contigs = m_graph->GetContigs();

	return 0;
}

// get two kmers from the line
void MapConverter::GetTwoKmer( string line, string &firstKmer, string &lastKmer ){
	vector<string> *temp = new vector<string>;
	Split( line, ",", temp );
	lastKmer = (*temp)[ 3 ];

	// get the first kmer, by removing the first space
	int spacePos = (*temp)[ 2 ].find( " " );
	firstKmer = (*temp)[ 2 ].substr( spacePos + 1, (*temp)[ 2 ].length() );
	//cout<<firstKmer<<"\t"<<lastKmer<<endl;

	temp->clear();
	delete temp;
}

// insert the kmer into map
void MapConverter::InsertKmerIntoMap( string kmer, map<string, list<int>* >* kmerList, int contigID, int ori ){
	// save head and tail
	map<string, list<int>* >::iterator firstIter = kmerList->find( kmer );
	if( firstIter == kmerList->end() ){
		// not exist
		list<int> *newList = new list<int>;
		newList->push_back( contigID );
		newList->push_back( ori );
		kmerList->insert( pair<string, list<int>* >(kmer, newList) );
	}
	else{
		firstIter->second->push_back( contigID );
		firstIter->second->push_back( ori );
	}
}

// find all overlap information of SOAPdenovo
void MapConverter::FindOverlapOfSOAP(){
	
	//cout<<"total first kmer: "<<m_headList->size()<<endl;
	//cout<<"total last kmer: "<<m_tailList->size()<<endl;

	int numOfEdge = 0;

	map<string, list<int>* >::iterator headIter = m_headList->begin();
	map<string, list<int>*>::iterator tailIter;
	for( headIter = m_headList->begin(); headIter != m_headList->end(); headIter++ ){
		string kmer = (*headIter).first;
		tailIter = m_tailList->find( (*headIter).first );
		if( tailIter != m_tailList->end() ){
			// there is overlap
			//cout<<"kmer is: "<<kmer<<endl;
			//cout<<"Tail: ";
			list<int>::iterator headListIter = (*headIter).second->begin();
			while( headListIter != (*headIter).second->end() ){
				Contig *headContig = m_graph->GetContigUsingPos( *headListIter );
				headListIter++;
				//int headOri = *headListIter;
				list<int>::iterator tailListIter = (*tailIter).second->begin();
				while( tailListIter != (*tailIter).second->end() ){
					Contig *tailContig = m_graph->GetContigUsingPos( *tailListIter );
					tailListIter++;
					//int tailOri = *tailListIter;
					if( tailContig->GetID() < headContig->GetID() ){
						// add edge
						numOfEdge++;

						// add to graph
						/*********/
				
					}
					
					tailListIter++;
				}
				
				headListIter++;
			}
		}
	}

	cout<<"There are "<<numOfEdge<<" edges\n";

	// delete all kmer map
	headIter = m_headList->begin();
	while( headIter != m_headList->end() ){
		(*headIter).second->clear();
		delete (*headIter).second;
		m_headList->erase( headIter++ );
	}

	tailIter = m_tailList->begin();
	while( tailIter != m_tailList->end() ){
		(*tailIter).second->clear();
		delete (*tailIter).second;
		m_tailList->erase( tailIter++ );
	}

	//PrintEdges();
}

// calculate mean and std using Interquartile Range
int MapConverter::IQR( vector<double> *dis, double &mean, double &std, string libName, double &percentageOfOutliers, double &numberOfUsedPairs,
		       PetLibrary *currentLibrary )
{
	ofstream disWriter( (Configure::OUTPUT_FOLDER + "distanceOfMatePairs_" + libName).c_str(), ios::out );
	if( disWriter == NULL ){
		cout<<"ERROR: Cannot open "<<(Configure::OUTPUT_FOLDER + "distanceOfMatePairs_" + libName)<<" file"<<endl;
		return -1;
	}

	cerr<<"distance size: "<<dis->size()<<endl;
	if( dis->size() == 0 ){
		cerr<<"No paired end reads are mapped onto the same contig\n";
		return 1;
	}
	
	// sort all the distances
	sort (dis->begin(), dis->end() );
	
	// find the value of 25% and 75% of the data
	int startID = (int)((double)dis->size() * 0.25);
	int endID = (int)((double)dis->size() * 0.75);

	// get the threshold for outliers
	double region = dis->at( endID ) - dis->at( startID );
	double lowerbound = dis->at( startID ) - 3 * region;
	double upperbound = dis->at( endID ) + 3 * region;
	
	// initialize the possibility of distance
	currentLibrary->InitPossibility( lowerbound, upperbound );

	// calculate mean and std excluding all outliers
	double sum = 0;
	double num = 0;
	string disStr = "";

	/*ofstream distributionWriter( (Configure::OUTPUT_FOLDER + "distribution").c_str() );
	if( distributionWriter == NULL ){
		cout<<"ERROR: Cannot open "<<(Configure::OUTPUT_FOLDER + "distribution")<<" file"<<endl;
		return -1;
	}
	string distributionResults = "distance\n";
	*/

	for( int i = 0; i < (int)dis->size(); i++ ){
		disStr += itos( dis->at( i ) ) + "\n";
		if( dis->at( i ) >= lowerbound && dis->at( i ) <= upperbound ){
			sum += dis->at( i );
			currentLibrary->AddOccuOfDistance( dis->at( i ) );
			//distributionResults += itos( dis->at( i ) ) + "\n";
			num++;
		}
	}
	
	currentLibrary->CalculatePossibility();

	//distributionWriter.write( distributionResults.c_str(), distributionResults.length() );
	//distributionWriter.close();

	disWriter.write( disStr.c_str(), disStr.length() );
	percentageOfOutliers = (1.0 - num/(double)dis->size()) * 100;
	numberOfUsedPairs = num;

	mean = sum / num;
	std = 0;
	// calculate std
	for( int i = 0; i < (int)dis->size(); i++ ){
		if( dis->at( i ) >= lowerbound && dis->at( i ) <= upperbound ){
			std += pow( dis->at( i ) - mean, 2 );
		}
	}

	std = std / num;
	std = pow( std, 0.5 );

	return 1;
}

// check if two lines represent a pair of reads
bool MapConverter::IsPair( string firstRead, string secondRead ){
	string firstName = firstRead.substr( 0, firstRead.length() -2 );
	string secondName = secondRead.substr( 0, secondRead.length() - 2 );

	if( firstName == secondName ){
		return true;
	}
	else
		return false;
}

// calculate the distance of two reads in the same contig
// including the read length
int MapConverter::CalculateReadDisOnSameContig( vector<string> *firstColumn, vector<string> *secondColumn ){
	int number[ 4 ];
	number[ 0 ] = atoi( firstColumn->at( 3 ).c_str() );
	number[ 1 ] = number[ 0 ] + firstColumn->at( 4 ).length();
	number[ 2 ] = atoi( secondColumn->at( 3 ).c_str() );
	number[ 3 ] = number[ 2 ] + secondColumn->at( 4 ).length();

	int min = number[ 0 ];
	int max = min;

	for( int i = 1; i < 4; i++ ){
		if( number[ i ] < min )
			min = number[ i ];
		if( number[ i ] > max )
			max = number[ i ];
	}

	return max - min;
}



// convert bowtie format into opera format
int MapConverter::ConvertBowtieFileMultiLib( PetLibrary *lib ){
	//int sumTime = 0;
	
	vector<string> *path = new vector<string>;
	Split( lib->GetFileName(), "/", path );
	vector<string> *names = new vector<string>;
	Split( path->at( path->size() - 1 ), ".", names );
	ofstream mapWriter( (Configure::OUTPUT_FOLDER + "pairedEndReads_" + names->at( 0 )).c_str(), ios::out );
	lib->SetNameWithoutPath( names->at( 0 ) );
	if( mapWriter == NULL ){
		cout<<"ERROR: Cannot open "<<(Configure::OUTPUT_FOLDER + "pairedEndReads_" + names->at( 0 ))<<" file"<<endl;
		names->clear();
		delete names;
		path->clear();
		delete path;
		return -1;
	}
	names->clear();
	delete names;
	path->clear();
	delete path;
	
	string head = "ID\tFirst Contig\tOrientation\tSecond Contig\tOrientation\tDistance\tStandard Deviation\n";
	mapWriter.write( head.c_str(), head.length() );

	list<mappingInfo*> *mapInfoList = lib->GetMappingInfo();
	
	// no pairing information
	if( mapInfoList->empty() )
		return 0;

	mappingInfo *preLine, *nextLine;
	preLine = mapInfoList->front();
	mapInfoList->pop_front();

	int petID = 1;			// record the paired-end reads ID

	//cerr<<"number of paired reads: "<<mapInfoList->size()<<endl;
	while( !mapInfoList->empty() ){
		nextLine = mapInfoList->front();
		mapInfoList->pop_front();
			
		string newLine = ConvertBowtieFormatMultiLib( preLine, nextLine, lib, petID++ ); 								     
		mapWriter.write( newLine.c_str(), newLine.size() );

		delete preLine;
		delete nextLine;

		if( !mapInfoList->empty() ){
			preLine = mapInfoList->front();
			mapInfoList->pop_front();
		}
		else{
			break;
		}
	}

	mapWriter.close();

	return 1;
}



// convert bowtie format
string MapConverter::ConvertBowtieFormatMultiLib( mappingInfo *firstRead, mappingInfo *secondRead, PetLibrary *currentLib, int id ){
	// get the current library
	string result = itos( id ) + "\t" + firstRead->contig->GetName() + "\t";

	Contig *contig1 = firstRead->contig;
	Contig *contig2 = secondRead->contig;
	int pos1 = m_graph->GetContigIndex( contig1->GetName() );
	int pos2 = m_graph->GetContigIndex( contig2->GetName() );
	double leftDis, rightDis;
	int dis;
	string firstOri, secondOri;

	if( currentLib->GetOri() == IN ){
		// orientation is ->...<- (+,-)
		if( firstRead->orientation == "+" ){
			leftDis = contig1->GetLength() - firstRead->position;
			firstOri = "+";
			if( secondRead->orientation == "+" ){
				// read ori is + +, contig ori should be + -
				result += "+\t" + secondRead->contig->GetName() + "\t-\t";
				secondOri = "-";
				rightDis = contig2->GetLength() - secondRead->position;
			}
			else{
				// read ori is + -; contig ori should be + +
				result += "+\t" + secondRead->contig->GetName() + "\t+\t";
				secondOri = "+";
				rightDis = secondRead->position + secondRead->readLength;
			}
		}
		else{
			leftDis = firstRead->position + firstRead->readLength;
			firstOri = "-";
			if( secondRead->orientation == "+" ){
				// read ori is - +; contig ori should be - -
				result += "-\t" + secondRead->contig->GetName() + "\t-\t";
				secondOri = "-";
				rightDis = contig2->GetLength() - secondRead->position;
			}
			else{
				// read ori is - -; contig ori should be - +
				result += "-\t" + secondRead->contig->GetName() + "\t+\t";
				secondOri = "+";
				rightDis = secondRead->position + secondRead->readLength;
			}
		}
	} // end dealing with orientation "in"
	else if( currentLib->GetOri() == OUT ){
		// orientation is <-...-> (-,+)
		if( firstRead->orientation == "+" ){
			leftDis = firstRead->position + firstRead->readLength;
			firstOri = "-";
			if( secondRead->orientation == "+" ){
				// read ori is + +, contig ori should be - +
				result += "-\t" + secondRead->contig->GetName() + "\t+\t";
				secondOri = "+";
				rightDis = secondRead->position + secondRead->readLength;
			}
			else{
				// read ori is + -; contig ori should be - -
				result += "-\t" + secondRead->contig->GetName() + "\t-\t";
				secondOri = "-";
				rightDis = contig2->GetLength() - secondRead->position;
			}
		}
		else{
			leftDis = contig1->GetLength() - firstRead->position;
			firstOri = "+";
			if( secondRead->orientation == "+" ){
				// read ori is - +; contig ori should be + +
				result += "+\t" + secondRead->contig->GetName() + "\t+\t";
				secondOri = "+";
				rightDis = secondRead->position + secondRead->readLength;
			}
			else{
				// read ori is - -; contig ori should be + -
				result += "+\t" + secondRead->contig->GetName() + "\t-\t";
				secondOri = "-";
				rightDis = contig2->GetLength() - secondRead->position;
			}
		}
	}
	else if( currentLib->GetOri() == FORWARD2 ){
		// orientation is ->...-> (+,+), 2nd read...1st read
		// swap read;
		mappingInfo *tempInfo = secondRead;
		secondRead = firstRead;
		firstRead = tempInfo;

		contig1 = firstRead->contig;
		contig2 = secondRead->contig;
		pos1 = m_graph->GetContigIndex( contig1->GetName() );
		pos2 = m_graph->GetContigIndex( contig2->GetName() );

		if( firstRead->orientation == "+" ){
			leftDis = contig1->GetLength() - firstRead->position;
			firstOri = "+";
			if( secondRead->orientation == "+" ){
				// read ori is + +, contig ori should be + +
				result += "+\t" + secondRead->contig->GetName() + "\t+\t";
				secondOri = "+";
				rightDis = secondRead->position + secondRead->readLength;
			}
			else{
				// read ori is + -; contig ori should be + -
				result += "+\t" + secondRead->contig->GetName() + "\t-\t";
				secondOri = "-";
				rightDis = contig2->GetLength() - secondRead->position;
			}
		}
		else{
			leftDis = firstRead->position + firstRead->readLength;
			firstOri = "-";
			if( secondRead->orientation == "+" ){
				// read ori is - +; contig ori should be - +
				result += "-\t" + secondRead->contig->GetName() + "\t+\t";
				secondOri = "+";
				rightDis = secondRead->position + secondRead->readLength;

			}
			else{
				// read ori is - -; contig ori should be - -
				result += "-\t" + secondRead->contig->GetName() + "\t-\t";
				secondOri = "-";
				rightDis = contig2->GetLength() - secondRead->position;
			}
		}
	}
	else if( currentLib->GetOri() == FORWARD1 ){
		// orientation is ->...-> (+,+), 1st read...2nd read
		if( firstRead->orientation == "+" ){
			leftDis = contig1->GetLength() - firstRead->position;
			firstOri = "+";
			if( secondRead->orientation == "+" ){
				// read ori is + +, contig ori should be + +
				result += "+\t" + secondRead->contig->GetName() + "\t+\t";
				secondOri = "+";
				rightDis = secondRead->position + secondRead->readLength;
			}
			else{
				// read ori is + -; contig ori should be + -
				result += "+\t" + secondRead->contig->GetName() + "\t-\t";
				secondOri = "-";
				rightDis = contig2->GetLength() - secondRead->position;
			}
		}
		else{
			leftDis = firstRead->position + firstRead->readLength;
			firstOri = "-";
			if( secondRead->orientation == "+" ){
				// read ori is - +; contig ori should be - +
				result += "-\t" + secondRead->contig->GetName() + "\t+\t";
				secondOri = "+";
				rightDis = secondRead->position + secondRead->readLength;

			}
			else{
				// read ori is - -; contig ori should be - -
				result += "-\t" + secondRead->contig->GetName() + "\t-\t";
				secondOri = "-";
				rightDis = contig2->GetLength() - secondRead->position;
			}
		}
	}
	else{
		leftDis = 1;
		rightDis = 1;
		cout<<"Dummy Case. Should not enter here";
	}
	// Any need to keep elseif and not make it else? ##?
	dis = (int) (currentLib->GetMean() - leftDis - rightDis);
	result += itos( dis ) + "\t" + itos( currentLib->GetStd() ) + "\n";
	if( pos1 > pos2 ){
		// change two contigs
		int tempPos = pos1;
		pos1 = pos2;
		pos2 = tempPos;
		string tempOri = firstOri;
		if( secondOri == "+" )
			firstOri = "-";
		else
			firstOri = "+";
		if( tempOri == "+" )
			secondOri = "-";
		else
			secondOri = "+";

	}
	// create a single pet object
	SinglePet *newPet = new SinglePet( id, pos1, firstOri, pos2, secondOri, dis, currentLib->GetStd() );
	//newPet->SetOriString( firstline + "\n" + secondline );
	//newPet->SetReadName( firstColumn->at( 0 ) );

	/*if( contig1->GetName() == "483457" && contig2->GetName() == "489913" ){
		cerr<<"find\n";
		cerr<<"distance : "<<dis<<endl;
		cerr<<newPet->GetOriString()<<endl;
	}
	*/

	// save to map
	int firstID = pos1;
	int secondID = pos2;
	if( firstOri == "-" )
		firstID = -firstID;
	if( secondOri == "-" )
		secondID = -secondID;

	// insert to pet map
	currentLib->InsertPair( firstID, secondID, newPet );

	return result;
}



bool less_singlePet( SinglePet*& p1, SinglePet*& p2 ){
	if( p1->GetStartID() < p2->GetStartID() ) 
		return true;
	else if( p1->GetStartID() > p2->GetStartID() )
		return false;
	else{
		if( p1->GetEndID() < p2->GetEndID() )
				return true;
			else if( p1->GetEndID() > p2->GetEndID() )
				return false;
		else{
			if( p1->GetStartOri() < p2->GetStartOri() )
				return true;
			else if( p1->GetStartOri() > p2->GetStartOri() )
				return false;
			else{
				if( p1->GetEndOri() < p2->GetEndOri() )
					return true;
				else if( p1->GetEndOri() > p2->GetEndOri() )
					return false;
				else
					return p1->GetDistance() < p2->GetDistance();
			}// finish checking second contig id
		} // finish checking first contig orientation
	}// finish checking first contig id
}


// start bundle
int MapConverter::Bundle(){
	ofstream clusterWriter( (Configure::OUTPUT_FOLDER + "clusters").c_str() );
	if( clusterWriter == NULL ){
		cout<<"ERROR: Cannot open "<<(Configure::OUTPUT_FOLDER + "clusters")<<" file"<<endl;
		return -1;
	}
	
	string head = "First Contig\tOrientation\tSecond Contig\tOrientation\tDistance\tStandard Deviation\tSize\n";
	clusterWriter.write( head.c_str(), head.length() );

	// traverse all group
	map<pair<int, int>, multiset<SinglePet*, less_distance>*>::iterator mapIter = m_singlePetsMap->begin();
	while( mapIter != m_singlePetsMap->end() ){
		// deal with current PET set
		string clusterString = BundleCluster( mapIter->second ); 
		clusterWriter.write( clusterString.c_str(), clusterString.length() );

		// proceed to next one
		delete mapIter->second;
		mapIter++;
	}

	m_singlePetsMap->clear();
	clusterWriter.close();
	return 1;
}

// create tables for gap correction
void MapConverter::CreateGapTables( PetLibrary *lib ){
	// clear the gap tables vector
	for( vector<GapCorrecter*>::iterator iter = m_gapTables.begin(); iter != m_gapTables.end(); iter++ ){
		delete *iter;
	}
	m_gapTables.clear();
	
	double maxDis = 0;
	maxDis = lib->GetMean() + Configure::STD_TIMES * lib->GetStd();
	if( lib->GetMaxDis() < maxDis )	
		maxDis = lib->GetMaxDis();

	//cerr<<"max dis: "<<lib->GetMaxDis()<<endl;
		
	for( int contigLengthSum = 500; contigLengthSum < maxDis + 500; contigLengthSum += 500 ){
		//cout<<"Constructing gap table for contig length: "<<contigLengthSum<<"...\n";

		GapCorrecter *gc = new GapCorrecter( maxDis, lib->GetMinDis(),  
						     contigLengthSum, 
						     lib->GetMean(), 
						     lib->GetStd(),
						     lib->GetReadLength(),
						     lib->GetPossibilities() );
		gc->CalculateTable();
		m_gapTables.push_back( gc );
	}
}

// start bundle
int MapConverter::BundleMultiLib( list<PetLibrary*> *libs, FILE *discardedEdgeFile, list<PetLibrary*> *superCluster, list<string> *conflictingEdges ){
	if( Configure::MAP_TYPE != OPERA ){
		for( list<PetLibrary*>::iterator libIter = libs->begin(); libIter != libs->end(); libIter++ ){
			PetLibrary *currentLib = *libIter;
			//Configure::CLUSTER_THRESHOLD = currentLib->GetThreshold();

			// create all the table for edge distance correction
			CreateGapTables( currentLib );
			//cerr<<"finish creating gap table\n";
			
			ofstream clusterWriter( (Configure::OUTPUT_FOLDER + "clusters_" + currentLib->GetNameWithoutPath() ).c_str() );
			ofstream clusterInfoWriter( (Configure::OUTPUT_FOLDER + "clustersInfo_" + currentLib->GetNameWithoutPath() ).c_str() );
			if( clusterWriter == NULL ){
				cout<<"ERROR: Cannot open "<<(Configure::OUTPUT_FOLDER + "clustersInfo_" + currentLib->GetNameWithoutPath() )<<" file"<<endl;
				return -1;
			}
			
			string head = "First Contig\tOrientation\tSecond Contig\tOrientation\tDistance\tStandard Deviation\tSize\n";
			clusterWriter.write( head.c_str(), head.length() );
			
			// traverse all group
			map<pair<int, int>, multiset<SinglePet*, lessDistance>*> *pairs = currentLib->GetSinglePetsMap();
			map<pair<int, int>, multiset<SinglePet*, lessDistance>*>::iterator mapIter = pairs->begin();
			while( mapIter != pairs->end() ){
				// deal with current PET set
				// correct the distance for each pet
				//cerr<<"Handling: "<<(*mapIter).first.first<<"\t"<<(*mapIter).first.second<<endl;
				
				// bundle those pets
				string clusterString = BundleClusterMultiLib( mapIter->second, currentLib, &clusterInfoWriter, discardedEdgeFile,
									      currentLib ); 
				//cerr<<"finish bundle one pair of contigs\n";
				clusterWriter.write( clusterString.c_str(), clusterString.length() );
				
				// proceed to next one
				//cerr<<"proceed to the next\n";
				delete mapIter->second;
				mapIter++;
				//cerr<<"finish one loop\n\n";
			}
			
			//cerr<<"finish one library\n";
			
			pairs->clear();
			clusterWriter.close();
			clusterInfoWriter.close();
			//cerr<<"ready to quit one loop\n\n";
			
			//cout<<"Total number of clusters: "<<currentLib->GetEdges()->size()<<endl;
		}
	}

	// bundle multiple libraries to create super clusters
	/*if( Configure::SIMULTANEOUSLY_HANDLE_LIBRARIES ){
		BundleSuperCluster( libs, superCluster, conflictingEdges );	
	}
	*/

	//cerr<<"start bundle multiple libraries\n";
	// bundle multiple libraries to create super clusters
	if( Configure::SIMULTANEOUSLY_HANDLE_LIBRARIES ){
		list<PetLibrary*> tempLibraries;

		// find all the small libraries and create one super cluster (<1kbp)
		list<PetLibrary*>::iterator libIter = libs->begin();
		while( libIter != libs->end() ){
			if( (*libIter)->GetMean() < Configure::SUPER_LIBRARY_THRESHOLD_1 ){
				tempLibraries.push_back( *libIter );	
				libIter = libs->erase( libIter );
			}
			else
				libIter++;
		}
		//cerr<<"bundle small libraries (<1kbp): "<<tempLibraries.size()<<endl;
		//cerr<<"small libraries: "<<tempLibraries.size()<<endl;
		if( tempLibraries.size() != 0 )
			BundleSuperCluster( &tempLibraries, superCluster, conflictingEdges, SMALL_LIB );

		// find all the medium libraries and create one super cluster (1kbp<= x < 10kbp)
		tempLibraries.clear();
		libIter = libs->begin();
		while( libIter != libs->end() ){
			if( (*libIter)->GetMean() >= Configure::SUPER_LIBRARY_THRESHOLD_1 
				&& (*libIter)->GetMean() < Configure::SUPER_LIBRARY_THRESHOLD_2 ){
				tempLibraries.push_back( *libIter );	
				libIter = libs->erase( libIter );
			}
			else
				libIter++;
		}
		//cerr<<"bundle medium libraries (1kbp to 10kbp): "<<tempLibraries.size()<<endl;
		//cerr<<“medium libraries: "<<tempLibraries.size()<<endl;
		if( tempLibraries.size() != 0 )
			BundleSuperCluster( &tempLibraries, superCluster, conflictingEdges, MEDIUM_LIB );

		// for all the other big libraries, create separate super cluster
		tempLibraries.clear();
		libIter = libs->begin();
		//cerr<<"start bundle big libraries: "<<libs->size()<<endl;
		// handle big libraries one by one
		while( libIter != libs->end() ){
			tempLibraries.push_back( *libIter );
			//cerr<<"bundle other library: "<<(*libIter)->GetFileName()<<endl;
			libIter = libs->erase( libIter );
			BundleSuperCluster( &tempLibraries, superCluster, conflictingEdges, BIG_LIB );
			//cerr<<"finish other library\n\n";
			tempLibraries.clear();
		}
		
		

		/*
		// handle all big libraries together
		while( libIter != libs->end() ){
			tempLibraries.push_back( *libIter );	
			libIter = libs->erase( libIter );
		}
		
		BundleSuperCluster( &tempLibraries, superCluster, conflictingEdges, false );
		tempLibraries.clear();
		*/
	}
	
	return 1;
}

// bundle multiple libraries to create super clusters
void MapConverter::BundleSuperCluster( list<PetLibrary*> *libs, list<PetLibrary*> *superCluster, list<string> *conflictingEdges, int libType ){
	m_ClustersMap->clear();

	int maxMean = 0;
	int maxStd = 0;
	int maxValue = 0;
	// combine all clusters together
	for( list<PetLibrary*>::iterator libIter = libs->begin(); libIter != libs->end(); libIter++ )
	{
		PetLibrary *currentLib = *libIter;

		if( currentLib->GetMean() + Configure::STD_TIMES * currentLib->GetStd() >= maxValue ){
			maxValue = currentLib->GetMean() + Configure::STD_TIMES * currentLib->GetStd();
			maxMean = currentLib->GetMean();
			maxStd = currentLib->GetStd();
		}

		list<PET*> *edges = currentLib->GetEdges();
		for( list<PET*>::iterator edgeIter = edges->begin(); edgeIter != edges->end(); edgeIter++ ){
			PET *currentEdge = (*edgeIter);
			string firstContig = currentEdge->GetStartContig()->GetName();
			if( currentEdge->GetStartContigOri() == MINUS )
				firstContig = "-" + firstContig;

			string secondContig = currentEdge->GetEndContig()->GetName();
			if( currentEdge->GetEndContigOri() == MINUS )
				secondContig = "-" + secondContig;
			
			pair<string, string> newPair( firstContig, secondContig );
			map<pair<string, string>, multiset<PET*, less_std>*>::iterator mapIter = m_ClustersMap->find( newPair );
			if( mapIter == m_ClustersMap->end() ){
				// new pair, create a new set
				multiset<PET*, less_std> *newSet = new multiset<PET*, less_std>;
				newSet->insert( currentEdge );
				m_ClustersMap->insert( pair<pair<string, string>, multiset<PET*, less_std>*>( newPair, newSet ) );
			}
			else{
				// insert the single pet to previous set
				mapIter->second->insert( currentEdge );
			}

		}
		
		// clear the list of edges of this library
		edges->clear();
	}
	
	// for each pair of contigs, 
	map<pair<string, string>, multiset<PET*, less_std>*>::iterator mapIter = m_ClustersMap->begin();
	string libraryName = "";
	if( libType == SMALL_LIB )
		libraryName = "super library (small libraries)";
	else if( libType == MEDIUM_LIB )
		libraryName = "super library (medium libraries)";
	else if( libs->size() > 1 )
		libraryName = "super library (big libraries)";
	else
		libraryName = (*(libs->begin()))->GetFileName();

	PetLibrary *newLibrary = new PetLibrary( libraryName );
	newLibrary->SetMean( maxMean );
	newLibrary->SetStd( maxStd );
	while( mapIter != m_ClustersMap->end() ){
		// deal with current PET set
		BundleOneSuperCluster( mapIter->second, newLibrary, conflictingEdges );
		/*if( mapIter->second->size() > 1 ){
			cout<<"smallest: "<<(*mapIter->second->begin())->GetStd()<<endl;
			cout<<"biggest: "<<(*mapIter->second->rbegin())->GetStd()<<endl<<endl;
			}*/
		
		// proceed to next one
		delete mapIter->second;
		mapIter++;
	}

	superCluster->push_back( newLibrary );

	m_ClustersMap->clear();

	for( list<PetLibrary*>::iterator libIter = libs->begin(); libIter != libs->end(); libIter++ )
	{
		delete (*libIter);
	}
	return;
}


// select the biggest cluster if there is a confliction
string MapConverter::BundleOneSuperCluster( multiset<PET*, less_std> *group, PetLibrary *lib, list<string> *conflictingEdges ){ 
	//PET *result;
	string finalCluster = "";

	int totalSize = 0;
	multiset<PET*, less_std>::iterator tempIter = group->begin();
	while( tempIter != group->end() ){
		totalSize += (*tempIter)->GetSize();
		tempIter++;
	}
	
	while( !group->empty() ){
		// find the edge with the maximum standard deviation
		// find all clusters within the range
		// do cluster
		// only accept the cluster containing more than half of the reads
		//int biggest = group->size() - 1;
		PET *biggestPet = (*group->rbegin());
		
		// initialize formula
		double p = 0;
		double q = 0;
		double newMean = 0;
		double newStd = 0;
		int size = 0;
		
		// calculate distance bound
		int lowerbound = biggestPet->GetDis() - Configure::STD_TIMES * biggestPet->GetStd();
		int upperbound = biggestPet->GetDis() + Configure::STD_TIMES * biggestPet->GetStd();
		string startName = biggestPet->GetStartContig()->GetName();
		//int startOri = biggestPet->GetStartContigOri();
		string endName = biggestPet->GetEndContig()->GetName();
		//int endOri = biggestPet->GetEndContigOri();
		bool startIsRepeat = biggestPet->GetStartContig()->IsRepeat();
		bool endIsRepeat = biggestPet->GetEndContig()->IsRepeat();
		//PET *superCluster = new PET();
		
		// get all edges in same bundle
		set<PET*> *subset = new set<PET*>;
		multiset<PET*, less_std>::iterator iter = group->begin();
		while( iter != group->end() ){
			//cerr<<(*iter)->GetOriString()<<endl;
			if( (*iter)->GetDis() > lowerbound && (*iter)->GetDis() < upperbound ){
				//cerr<<"in range\n";
				size += (*iter)->GetSize();
				subset->insert( *iter );
				multiset<PET*, less_std>::iterator preIter = iter;
				iter++;
				group->erase( preIter );
			}
			else
				iter++;
		}
		
		// for edges between unique contigs, only keep the biggest bundle
		// also discard all multiple edges between uniqe contigs for polyploid sequence
		// also discard bundles smaller than threshold
		if( (!startIsRepeat && !endIsRepeat && size <= (double) totalSize / 2 ) 
			|| (!startIsRepeat && !endIsRepeat && Configure::PLOIDY > 1 && size != totalSize )
			|| size < Configure::CLUSTER_THRESHOLD ){
			// ignore current subset
			//cout<<"Warning: There is a small super cluster between a pair of contigs: "<<biggestPet->GetStartContig()->GetName()
			//    <<"\t"<<biggestPet->GetEndContig()->GetName()<<endl;
			//cout<<"size: "<<size<<" total size: "<<totalSize<<endl;
			set<PET*>::iterator  iter = subset->begin();
			while( iter != subset->end() ){
				conflictingEdges->push_back( (*iter)->GetOriString() );
				delete *iter;
				iter++;
			}
		}
		else{
			// form one super cluster
			size = 0;
			PET *superCluster = new PET();
			//cout<<"There is only one super cluster\n";
			set<PET*>::iterator iter = subset->begin();
			Contig* startContig = (*iter)->GetStartContig();
			Contig* endContig = (*iter)->GetEndContig();
			int startOri = (*iter)->GetStartContigOri();
			int endOri = (*iter)->GetEndContigOri();
			while( iter != subset->end() ){
				// calculate the value
				//p += (*iter)->GetDis() / pow( (double)(*iter)->GetStd(), 2 );
				//q += 1 / pow( (double)(*iter)->GetStd(), 2 );
				
				// new way to calculate the values
				p += (double) (*iter)->GetDis() * (double)(*iter)->GetSize() / pow( pow( (*iter)->GetSize(), 0.5 ) * (*iter)->GetStd(), 2 );
				q += (double) (*iter)->GetSize() / pow( pow( (*iter)->GetSize(), 0.5 ) * (*iter)->GetStd(), 2 );
				
				// erase from the group
				superCluster->AddOneClusterToSuperCluster( *iter );
				size += (*iter)->GetSize();
				iter++;
			}

			// calculate new mean and std
			newMean = p / q;
			newStd = 1.0 / pow( q, 0.5 );
			// consider overlap between contigs
			if( (int)newStd == 0 )
				newStd = 1;
			
			superCluster->SetProperties( startContig, startOri, 
					     endContig, endOri, 
					     (int)(newMean+0.5), (int)(newStd+0.5), size );


			lib->AddCluster( superCluster );
		
			// add to multiple libraries
			m_graph->AddEdgeMultiLib( superCluster );
		}

		subset->clear();
		delete subset;
	}
	return finalCluster;
}


// bundle a certain cluster
// return all the cluster string
string MapConverter::BundleCluster( multiset<SinglePet*, less_distance> *group ){
	PET *result;
	string finalCluster = "";

	if( group->empty() )
		return "";

	while( !group->empty() ){
		// find the median single pet
		int median = group->size() / 2 + 1;

		multiset<SinglePet*, less_distance>::iterator iter = group->begin();
		int counter = 1;
		while( counter < median ){
			iter++;
			counter++;
		}
		SinglePet *medianPet = *iter;

		// initialize formula
		double p = 0;
		double q = 0;
		double newMean = 0;
		double newStd = 0;
		int size = 0;
		
		// calculate distance bound
		int lowerbound = medianPet->GetDistance() - Configure::STD_TIMES * medianPet->GetStd();
		int upperbound = medianPet->GetDistance() + Configure::STD_TIMES * medianPet->GetStd();
		int startID = medianPet->GetStartID();
		int startOri = medianPet->GetStartOri();
		int endID = medianPet->GetEndID();
		int endOri = medianPet->GetEndOri();

		// get all edges in same bundle
		iter = group->begin();
		while( iter != group->end() && (*iter)->GetDistance() < upperbound ){
			if( (*iter)->GetDistance() > lowerbound && (*iter)->GetDistance() < upperbound ){
				// in the bundle
				size++;

				// calculate the value
				p += (*iter)->GetDistance() / pow( (double)(*iter)->GetStd(), 2 );
				q += 1 / pow( (double)(*iter)->GetStd(), 2 );

				// remove that list
				delete (*iter);
				group->erase( iter++ );
			}
			else{
				iter++;
			}
		}

		// calculate new mean and std
		newMean = p / q;
		newStd = 1 / pow( q, 0.5 );
		// consider overlap between contigs
		if( newMean < 0 )
			newMean += (-Configure::MIN_DIS);

		if( (int)newStd == 0 )
			newStd = 1;
		

		// add to graph
		result = new PET( m_graph->GetContig( startID ), startOri, m_graph->GetContig( endID ), endOri, (int) newMean, (int) newStd, size );
		finalCluster += result->ToString() + "\n";
		if( size >= Configure::CLUSTER_THRESHOLD
		    && newMean + Configure::STD_TIMES * newStd >= 0
		    && newMean - Configure::STD_TIMES * newStd <= Configure::LIB_MEAN + Configure::LIB_STD * Configure::STD_TIMES )
			m_graph->AddEdge( result );
		else
			delete result;
	}

	return finalCluster;
}

// correct the distance of clusters
/*
void MapConverter::CorrectDistance( multiset<SinglePet*, lessDistance> *group, PetLibrary *lib ){
	if( group->empty() )
		return;

	if( m_graph->GetContig( (*group->begin())->GetStartID() ) ); // "If" discards return value
	if( m_graph->GetContig( (*group->begin())->GetEndID() ) ); // "If" discards return value

	// create the table of g and g_head
	//int libUpperbound = Configure::UPPERBOUND; 
	//int contigSumLength = startContig->GetLength() + endContig->GetLength();

	

	// for each pet, find the closest value of g_head, and change it to the corresponding g

}
*/

// bundle a certain cluster
string MapConverter::BundleClusterMultiLib( multiset<SinglePet*, lessDistance> *group, PetLibrary *lib, ofstream *clusterInfoFile, FILE *discardedEdgeFile,
					    PetLibrary *currentLibrary ){
	PET *result;
	string finalCluster = "";

	if( group->empty() )
		return "";

	//bool isSpecialEdge = false;

	//cerr<<"start bundle\n";
	string clusterInfo = "";
	while( !group->empty() ){
		// find the median single pet
		//cerr<<"try to find median\n";
		int median = group->size() / 2 + 1;
		bool isEdge = false;	

		multiset<SinglePet*, lessDistance>::iterator iter = group->begin();
		int counter = 1;
		while( counter < median ){
			iter++;
			counter++;
		}
		
		SinglePet *medianPet = *iter;

		// initialize formula
		double p = 0;
		double q = 0;
		double newMean = 0;
		double newStd = 0;
		int size = 0;
		
		// calculate distance bound
		int lowerbound = medianPet->GetDistance() - Configure::STD_TIMES * medianPet->GetStd();
		int upperbound = medianPet->GetDistance() + Configure::STD_TIMES * medianPet->GetStd();
		int startID = medianPet->GetStartID();
		int startOri = medianPet->GetStartOri();
		int endID = medianPet->GetEndID();
		int endOri = medianPet->GetEndOri();

		//cerr<<"lower bound: "<<lowerbound<<"\tupper bound: "<<upperbound<<endl;
		//cerr<<"median reads distance: "<<medianPet->GetDistance()<<endl;

		// get all edges in same bundle
		iter = group->begin();
#ifdef ANCHOR
		list<SinglePet*> *bundlePet = new list<SinglePet*>;
#endif

		//double sum = 0;
		//double num = 0;
		while( iter != group->end() && (*iter)->GetDistance() <= upperbound ){
			if( (*iter)->GetDistance() >= lowerbound && (*iter)->GetDistance() <= upperbound ){
				clusterInfo += (*iter)->GetReadName().substr( 0, (*iter)->GetReadName().length() - 2 ) + "\t";
				// in the bundle
				size++;

				// calculate the value
				p += (double)(*iter)->GetDistance() / pow( (double)(*iter)->GetStd(), 2 );
				q += 1.0 / pow( (double)(*iter)->GetStd(), 2 );
				//sum+= (*iter)->GetDistance();
				//num++;

				// remove that list
#ifdef ANCHOR
				bundlePet->push_back( *iter );
#endif

#ifndef ANCHOR
				delete (*iter);
#endif
				group->erase( iter++ );
			}
			else{
				iter++;
			}
		}
		
		if( isEdge ){
			cerr<<"size is "<<size<<endl;
		}
		
#ifdef ANCHOR
		string anchorString = CalculateAnchorRegion( bundlePet );
		bundlePet->clear();
		delete bundlePet;
#endif

		// calculate new mean and std
		newMean = p / q;
		newStd = 1.0 / pow( q, 0.5 );
		// consider overlap between contigs
		//if( newMean < 0 )
		//	newMean += (-Configure::MIN_DIS);

		if( (int)(newStd+0.5) == 0 )
			newStd = 1;

		// add to graph
		/*********/
		result = new PET( m_graph->GetContig( startID ), startOri, m_graph->GetContig( endID ), endOri, (int) (newMean+0.5), (int) (newStd+0.5), size );

		// refine the mean of this cluster
		//cerr<<"start refining mean\n";
		//cerr<<"old mean is "<<newMean<<endl;
		//cerr<<"two contigs are: "<<result->GetStartContig()->GetName()<<"\t"<<result->GetEndContig()->GetName()<<endl;
		RefineMean( result, result->GetStartContig(), result->GetEndContig() );
		//cerr<<"end refining mean\n\n\n";

#ifdef ANCHOR
		result->AddAnchorString( anchorString );
		//cout<<result->GetOriString()<<endl;
		finalCluster += result->ToString() + anchorString + "\n";
#endif

#ifndef ANCHOR
		finalCluster += result->ToString() + "\n";
#endif
		clusterInfo += "\n";
		clusterInfo = result->ToString() + "\n" + clusterInfo;
		clusterInfoFile->write( clusterInfo.c_str(), clusterInfo.length() );
		clusterInfo = "";

		int startCopyNumber = (int)(result->GetStartContig()->GetCov()/Configure::HAPLOID_COVERAGE + 0.5);
		int endCopyNumber = (int)(result->GetEndContig()->GetCov()/Configure::HAPLOID_COVERAGE + 0.5);

		if( size >= Configure::CLUSTER_THRESHOLD &&
			newMean + Configure::STD_TIMES * newStd + Configure::KMER - 1 >= 0
			&& newMean - Configure::STD_TIMES * newStd <= lib->GetMean() + lib->GetStd() * Configure::STD_TIMES 
		    && (result->GetStartContig()->GetContigType() != SMALL 
		    	&& result->GetEndContig()->GetContigType() != SMALL ) ){

			if( Configure::FILTER_REPEAT ){
				if( result->GetStartContig()->IsRepeat() || result->GetEndContig()->IsRepeat() ){
					// print to discarded edge file
					//cerr<<"start delete edge\n";
					fprintf( discardedEdgeFile, "%s\n", result->ToString().c_str() );
					delete result;
					//cerr<<"finish discard edge\n";
					continue;
				}
			}
			else{
				// remove the edges between repeats
				if( result->GetStartContig()->IsRepeat() && result->GetEndContig()->IsRepeat() ){
					// print to discarded edge file
					//cerr<<"start delete edge\n";
					fprintf( discardedEdgeFile, "%s\n", result->ToString().c_str() );
					delete result;
					//cerr<<"finish discard edge\n";
					continue;
				}
			}

			// if neither of contigs is repeat, check if they have the same copy number, 
			// and the difference of coverage is small
			if( !result->GetStartContig()->IsRepeat() && !result->GetEndContig()->IsRepeat() ){
				if( Configure::PLOIDY > 2 && ( startCopyNumber != endCopyNumber 
					|| abs( result->GetStartContig()->GetCov() - result->GetEndContig()->GetCov() ) 
					   > 0.5 * Configure::HAPLOID_COVERAGE ) ){
					// print to discarded edge file
					//cerr<<"start delete edge\n";
					fprintf( discardedEdgeFile, "%s\n", result->ToString().c_str() );
					delete result;
					//cerr<<"finish discard edge\n";
					continue;
				}
			}
			
			//cerr<<"start add edge to multi library\n";
			if( !Configure::SIMULTANEOUSLY_HANDLE_LIBRARIES )
				m_graph->AddEdgeMultiLib( result );
			//cerr<<"finish add edge to multi library\n";
			
			//cerr<<"start add cluster\n";
			lib->AddCluster( result );
			//cerr<<"finish add cluster\n";
			
			//cout<<result->ToString()<<endl;
		}
		else{
			// print to discared edge file
			//cerr<<"directly delete edge\n";
			fprintf( discardedEdgeFile, "%s\n", result->ToString().c_str() );
			delete result;
			//cerr<<"end diretly delete edge\n";
		}
	}

	//cerr<<"Prepare to return\n";
	return finalCluster;
}

// refine the mean of a cluster
void MapConverter::RefineMean( PET *e, Contig *c1, Contig *c2 ){
	// find the sum length of contigs
	int sumLength = c1->GetLength() + c2->GetLength();
	
	// find the corresponding table
	int id = int( sumLength / 500 - 1 );
	if( id < 0 )
		id = 0;

	if( id >= (int)m_gapTables.size() )
		id = m_gapTables.size() - 1;

	//cerr<<"total table is "<<m_gapTables.size()<<endl;
	//cerr<<"find "<<id<<" table\n";
  	
	GapCorrecter *currentGC = m_gapTables.at( id );

	//cerr<<"distance is "<<e->GetDis()<<endl;
   
	// find the real gap size
	// modify the distance of each paired-end reads
	e->SetDis( currentGC->GetRealGapSize( e->GetDis() ) );
}

// debugging purpose
string MapConverter::CalculateAnchorRegion( list<SinglePet*> *pets ){
	string results = "";
	for( int i = 0; i < 2; i++ ){
		// get the contig name
		//vector<string> *contents = new vector<string>;
		string contigName;
		if( i == 0 )
			contigName = m_graph->GetContig( pets->front()->GetStartID() )->GetName();
		else
			contigName = m_graph->GetContig( pets->front()->GetEndID() )->GetName();

		list<int> *regions = new list<int>;

		// for each pair of mapping, select the corresponding contig
		for( list<SinglePet*>::iterator iter = pets->begin(); iter != pets->end(); iter++ ){
			vector<string> *lines = new vector<string>;
			Split( (*iter)->GetOriString(), "\n", lines );
			int startPos = 0; int endPos = 0;
			for( int lineID = 0; lineID < 2; lineID++ ){
				vector<string> *contents = new vector<string>;
				Split( lines->at( lineID ), "\t", contents );
				if( contents->at( 2 ) == contigName ){
					// find the position
					//cout<<lines->at( lineID );
					if( Configure::MAP_TYPE == BOWTIE ){
						startPos = atoi( contents->at( 3 ).c_str() );
						endPos = startPos + contents->at( 4 ).length() - 1;
					}
					else if( Configure::MAP_TYPE == SAM ){
						startPos = atoi( contents->at( 3 ).c_str() );
						endPos = startPos + contents->at( 9 ).length() - 1;
					}
					//cout<<endl<<startPos<<"\t"<<endPos<<endl;
				}
				contents->clear();
				delete contents;
			}

			// save anchor
			list<int>::iterator regionIter = regions->begin();
			bool insert = false;
			while( regionIter != regions->end() ){
				list<int>::iterator start = regionIter;
				list<int>::iterator end = regionIter;
				end++;
				if( startPos <= *start ){
					regions->insert( regionIter, startPos );
					regions->insert( regionIter, endPos );
					insert = true;
					break;
				}
				regionIter++;
				regionIter++;
			}
			if( !insert ){
				regions->push_back( startPos );
				regions->push_back( endPos );
			}

			lines->clear();
			delete lines;
		}

		// update anchor
		//	if( pets->size() > 50 ){
		//	cout<<contigName<<endl;
		//cout<<pets->size()<<endl;
		//cout<<regions->size()<<endl;
		//}
		
		list<int>::iterator regionIter = regions->begin();
		int anchorLength = 0;
		int num = 0;
		while( regionIter != regions->end() ){	
			list<int>::iterator start1 = regionIter;
			list<int>::iterator end1 = start1;  end1++;
			list<int>::iterator start2 = end1; start2++;
			if( start2 == regions->end() ){
				anchorLength += *end1 - *start1 + 1;
				break;
			}
			list<int>::iterator end2 = start2; end2++;

			// check if the first pair overlaps with the second pair
			if( *start2 <= *end1 ){
				// if yes, combine those two pairs
				if( *end1 < *end2 )
					*end1 = *end2;

				regions->erase( start2 );
				regions->erase( end2 );
			}
			else{
				// if not, move the iterator
				regionIter++;
				regionIter++;
				anchorLength += *end1 - *start1 + 1;
				num++;
			}
		}
		
		//	if( pets->size() > 50 ){
			//	cout<<"anchor region length: "<<anchorLength<<endl<<endl;
		//}

		if( anchorLength == 0 ){
			cout<<"Find a 0 region: Total region is "<<num<<endl;
		}
		results += "\t" + itos( anchorLength );
		
		regions->clear();
		delete regions;
	}

	for( list<SinglePet*>::iterator iter = pets->begin(); iter != pets->end(); iter++ ){
		delete *iter;
	}

	return results;
}

// analyze opera edge file
int MapConverter::AnalyzeOpera( string fileName ){
	ifstream edgeReader( fileName.c_str() );

	if( edgeReader == NULL )
	{
		cout<<"ERROR: Cannot open "<<fileName<<" file"<<endl;
		return -1;
	}

	string line;
	vector<string> *currentPet = new vector<string>;
	while( getline( edgeReader, line ) ){
		Split( line, "\t", currentPet );

		if( currentPet->size() == 7 ){
			// add to graph
			Contig *firstContig = m_graph->GetContig( m_graph->GetContigIndex( (*currentPet)[ 0 ] ) );
			Contig *secondContig = m_graph->GetContig( m_graph->GetContigIndex( (*currentPet)[ 2 ] ) );
			int firstOri = PLUS;
			if( (*currentPet)[ 1 ] == "-" )
				firstOri = MINUS;
			int secondOri = PLUS;
			if( (*currentPet)[ 3 ] == "-" )
				secondOri = MINUS;

			PET *result = new PET( firstContig, firstOri, secondContig, secondOri,
				atoi( (*currentPet)[ 4 ].c_str() ), atoi( (*currentPet)[ 5 ].c_str() ),
				atoi( (*currentPet)[ 6 ].c_str() ) );
			if( atoi( (*currentPet)[ 6 ].c_str() ) >= Configure::CLUSTER_THRESHOLD ){
				if( result->GetDis() + Configure::STD_TIMES * result->GetStd() > 0  
					&& result->GetDis() - Configure::STD_TIMES <= Configure::LIB_MEAN + Configure::LIB_STD * Configure::STD_TIMES )
					m_graph->AddEdge( result );
			}
			else{
				delete result;
			}
		}
		currentPet->clear();
	}
	delete currentPet;
	return 1;
}

// analyze opera edge file
int MapConverter::AnalyzeOperaMultiLib( list<PetLibrary*> *libs ){
	// get all cluster files and handle them one by one
	//cerr<<"Total file in AnalyzeOperaMultiLib is "<<Configure::MULTI_LIB_INFO->size()<<endl;
	for( int i = 0; i < (int)Configure::MULTI_LIB_INFO->size(); i++ )
	{
		if(Configure::MULTI_LIB_INFO->at( i )-> GetMapType() == OPERA){
			cerr<<"Analyzing file: "<<Configure::MULTI_LIB_INFO->at( i )->GetFileName()<<endl;
			
			//cout<<i<<": threshold: "<<Configure::MULTI_LIB_INFO->at( i )->GetThreshold()<<endl;
			ifstream edgeReader( Configure::MULTI_LIB_INFO->at( i )->GetFileName().c_str() );

			if( edgeReader == NULL )
				{
					cout<<"ERROR: Cannot open "<<Configure::MULTI_LIB_INFO->at( i )->GetFileName()<<" file"<<endl;
					return -1;
				}
		
			PetLibrary *lib = new PetLibrary( Configure::MULTI_LIB_INFO->at( i )->GetFileName() );
			libs->push_back( lib );
			lib->SetMean( Configure::MULTI_LIB_INFO->at( i )->GetMean() );
			lib->SetStd( Configure::MULTI_LIB_INFO->at( i )->GetStd() );
			//lib->SetThreshold( Configure::MULTI_LIB_INFO->at( i )->GetThreshold() );
			//Configure::CLUSTER_THRESHOLD = lib->GetThreshold();

			string line;
			vector<string> *currentPet = new vector<string>;
			getline( edgeReader, line );
			while( getline( edgeReader, line ) ){
				Split( line, "\t", currentPet );

				//if( (*currentPet)[ 0 ] == "10639727" || (*currentPet)[ 2 ] == "10639727" ){
				//	cerr<<"original edge: "<<line<<endl;
				//}

				if( atoi( (*currentPet)[ 6 ].c_str() ) < Configure::CLUSTER_THRESHOLD )
					continue;

				// add to graph
				int firstID = m_graph->GetContigIndex( (*currentPet)[ 0 ] );
				int secondID = m_graph->GetContigIndex( (*currentPet)[ 2 ] );

				// if one of contigs does not exist
				if( firstID == -1 || secondID == -1 )
					continue;

				Contig *firstContig = m_graph->GetContig( firstID );
				Contig *secondContig = m_graph->GetContig( secondID );

				//if( firstContig->GetName() == "10639727" || secondContig->GetName() == "10639727" ){
				//	cerr<<"    edge with large threshold: "<<line<<endl;
				//}

				int firstOri = PLUS;
				if( (*currentPet)[ 1 ] == "-" )
					firstOri = MINUS;
				int secondOri = PLUS;
				if( (*currentPet)[ 3 ] == "-" )
					secondOri = MINUS;
			
				int edgeStd = atoi( (*currentPet)[ 5 ].c_str() );
				if( edgeStd == 0 )
					edgeStd = 1;

				PET *result = new PET( firstContig, firstOri, secondContig, secondOri,
						       atoi( (*currentPet)[ 4 ].c_str() ), edgeStd,
						       atoi( (*currentPet)[ 6 ].c_str() ) );
			
				if( result->GetSize() >= Configure::CLUSTER_THRESHOLD && 
				    result->GetDis() + Configure::STD_TIMES * result->GetStd()  + Configure::KMER - 1 >= 0  
				    && result->GetDis() - Configure::STD_TIMES * result->GetStd() <= lib->GetMean() + lib->GetStd() * Configure::STD_TIMES 
				    && (result->GetStartContig()->GetContigType() != SMALL && result->GetEndContig()->GetContigType() != SMALL ) ){
					if( Configure::FILTER_REPEAT ){
						if( result->GetStartContig()->IsRepeat() || result->GetEndContig()->IsRepeat() ){
							// print to discared edge file
							//cerr<<"filter repeat edges: "<<line<<endl;
							delete result;
							continue;
						}
					}
					else{
						// remove the edges between repeats
						if( result->GetStartContig()->IsRepeat() && result->GetEndContig()->IsRepeat() ){
							// print to discared edge file
							//cerr<<"filter edges between repeats: "<<line<<endl;;
							delete result;
							continue;
						}
					}
			       

					if( !Configure::SIMULTANEOUSLY_HANDLE_LIBRARIES )
						m_graph->AddEdgeMultiLib( result );
				
					lib->AddCluster( result );
					//cerr<<"add edge\t"<<result->ToString()<<endl;
					//if( firstContig->GetName() == "10639727" || secondContig->GetName() == "10639727" ){
					//	cerr<<"    save edge in file: "<<line<<endl;
					//}

				
				}
				else{
					//cerr<<"edge is not correct: "<<line<<endl;
					delete result;
				}
				currentPet->clear();
			}
			delete currentPet;
		}
	}

	// sort all libraries according to library means
	libs->sort( LibSort );

	return 1;
}
