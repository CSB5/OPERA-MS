#include "opera.h"
// optimal scaffolding

opera::opera(void)
{
	m_graph = new Graph();
	m_activeRegion = new list<Contig*>;		
	m_happyDanglingEdges = new list<PET*>;
	m_unhappyDanglingEdges = new list<PET*>;
	m_unassignedNodes = new list< pair<Contig*, int> >;
	m_unhappyEdgeString = new list<string>;
	m_visitedTree = new Tree();

	m_conflictingEdges = new list<string>;

	// for multi libraries
	m_libraries = new list<PetLibrary*>;
	m_superLibrary = new list<PetLibrary*>;

	m_graphIDIncreaseThreshold = 0;
	m_realMapping = new map<string, reference_mapping*>;

	m_totalPartialScaffold = 0;
        m_numOfTopUnassignedNodes = 0;
	m_repeatTraceBack = NULL;
}

// initialize for another run
void opera::Initialize(){
	delete m_visitedTree;
	m_visitedTree = new Tree();
	m_activeRegion->clear();;		
	m_happyDanglingEdges->clear();
	m_unhappyDanglingEdges->clear();
	m_unassignedNodes->clear();
}

opera::~opera(void)
{
	delete m_graph;
	m_activeRegion->clear();
	delete m_activeRegion;
	m_happyDanglingEdges->clear();
	delete m_happyDanglingEdges;
	m_unhappyDanglingEdges->clear();
	delete m_unhappyDanglingEdges;
	m_unassignedNodes->clear();
	delete m_unassignedNodes;

	delete m_unhappyEdgeString;
	delete m_visitedTree;
	fclose( logFile );
	fclose( m_discardEdgeFile );

	for( list<PetLibrary*>::iterator iter = m_libraries->begin(); iter != m_libraries->end(); iter++ ){
		delete *iter;
	}
	m_libraries->clear();
	delete m_libraries;
	for( list<PetLibrary*>::iterator iter = m_superLibrary->begin(); iter != m_superLibrary->end(); iter++ ){
		delete *iter;
	}
	m_superLibrary->clear();
	delete m_superLibrary;

	for( map<string, reference_mapping*>::iterator iter= m_realMapping->begin(); iter != m_realMapping->end(); iter++ ){
		delete (*iter).second;
	}
	m_realMapping->clear();
	delete m_realMapping;

	m_conflictingEdges->clear();
	delete m_conflictingEdges;
}

// clear scaffold variables
void opera::Clear(){
	m_activeRegion->clear();		
	m_happyDanglingEdges->clear();
	m_unhappyDanglingEdges->clear();
	m_unassignedNodes->clear();
	m_visitedTree->Clear();
}

// check if the names have confliction
// return true if has confliction
// return false if not
bool opera::CheckNameConfliction(){
	ifstream contigReader( Configure::CONTIG_FILE.c_str() );

	if( contigReader == NULL )
	{
		cout<<"ERROR: Cannot open "<<Configure::CONTIG_FILE<<" file"<<endl;
		return true;
	}
		
	string line;					
			
	getline( contigReader, line );
	contigReader.close();
	
	// check the contig name
	if( line.substr( 1, Configure::SCAFFOLD_PREFEX.length() ) == Configure::SCAFFOLD_PREFEX )
		return true;
	else
		return false;
}

// check if all files exist
bool opera::CheckFileExist(){
	// check contig file
	ifstream contigReader( Configure::CONTIG_FILE.c_str() );

	if( contigReader == NULL )
	{
		cout<<"ERROR: Cannot open "<<Configure::CONTIG_FILE<<" file"<<endl;
		return false;
	}
	
	// check mapping file
	for( int i = 0; i < (int) Configure::MULTI_LIB_INFO->size(); i++ )
	{
		ifstream mapReader( Configure::MULTI_LIB_INFO->at( i )->GetFileName().c_str() );

		if( mapReader == NULL )
		{
			cout<<"ERROR: Cannot open "<<Configure::MULTI_LIB_INFO->at( i )->GetFileName()<<" file"<<endl;
			return false;
		}

		mapReader.close();
	}
	contigReader.close();
	return true;


}

// read real position file
int opera::ReadRealPositionFile( string fileName ){
	ifstream realPositionReader( fileName.c_str() );

	//cerr<<"analyzing position file\n";

	if( realPositionReader == NULL )	{
		cout<<"ERROR: Cannot open "<<fileName<<" file"<<endl;
		return -1;
	}

/*
	string line;					
	vector<string> *contents = new vector<string>;		
	
	int id = 0;
	string previousReferenceName = "";
	while( getline( realPositionReader, line ) ){
		Split( line, "\t", contents );
		if( contents->at( 3 ) != previousReferenceName ){
			id = 0;
			previousReferenceName = contents->at( 3 );
		}

		int index = m_graph->GetContigIndex( contents->at( 0 ) );
		if( index != -1 ){
			Contig *c = m_graph->GetContigUsingPos( index );
			c->SetReferenceName( contents->at( 3 ) );
			c->SetReferenceStartPos( atoi( contents->at( 4 ).c_str() ) );
			c->SetReferenceEndPos( atoi( contents->at( 5 ).c_str() ) );
			c->SetReferenceIndex( id++ );
			if( contents->at( 2 ) == "BE" )
				c->SetReferenceOri( PLUS );
			else
				c->SetReferenceOri( MINUS );
		}

		// Old commented part
		reference_mapping *newMapping = new reference_mapping;
		newMapping->m_contigName = contents->at( 0 );
		newMapping->m_referenceGenomeName = contents->at( 3 );
		newMapping->m_referenceStartPos = atoi( contents->at( 4 ).c_str() ); 
		newMapping->m_referenceEndPos = atoi( contents->at( 5 ).c_str() ); 
		if( contents->at( 2 ) == "BE" )
			newMapping->m_referenceOri = PLUS;
		else
			newMapping->m_referenceOri = MINUS;
		newMapping->m_referenceIndex = id++;

		m_realMapping->insert( pair<string, reference_mapping*> ( contents->at( 0 ), newMapping ) );
	}

	contents->clear();
	delete contents;
*/	
	realPositionReader.close();

	return 1;

}

// comparison, not case sensitive.
bool compare_real_position( Contig *first, Contig *second){
	return first->GetReferenceStartPos() < second->GetReferenceStartPos();
}

// label the order of contigs in the subgraph according to the real mapping
void opera::LabelContigsUsingMapping( list<Contig*> *subgraph ){
	// sort the contigs according to the mapping position
	subgraph->sort( compare_real_position );

	// label the contigs
	int id = 0;
	for( list<Contig*>::iterator iter = subgraph->begin(); iter != subgraph->end(); iter++ ){
		(*iter)->SetSubgraphReferenceIndex( id++ );
	}
}

// start doing scaffolding
// return -1 if failed
int opera::StartScaffold(){

	#ifdef DEBUG
		int finishedScaffold = 0;
	#endif

	//int abandendGraph = 0;
	//int abandendNode = 0;
	bool ifIncreaseThreshold;
	int numOfSpecial = 1;
	m_totalNumberOfNodes = m_graph->GetNumOfContigs();

	while( m_graph->HasEdges() ){
		// find subgraph
		int numOfContig = 0;
		int numOfBorderContig = 0;
		int minClusterSize = 99999;			// the minimum cluster size in this subgraph
		list<Contig*> *subgraph = m_graph->FindSubgraph( numOfContig, numOfBorderContig, minClusterSize );
		m_subgraphSize = subgraph->size();
		m_totalNumberOfSubgraph++;

		bool printGraph = false;
		for( list<Contig*>::iterator iter = subgraph->begin(); iter != subgraph->end(); iter++ ){
			if( (*iter)->GetName() == "10639727" ){
				printGraph = true;
			}
		}
		//cout<<"minClusterSize is "<<minClusterSize<<endl;
		//cout<<"number of contigs is "<<numOfContig<<endl;

		// set all contigs unassigned nodes iterator to null
		//cerr<<"Scaffolding a subgraph with "<<subgraph->size()<<" nodes\n";
#ifdef LOG
		fprintf( logFile, "Scaffolding a subgraph with %d nodes...\n", (int) subgraph->size() );
		fflush( logFile );
		//cout<<"Scaffolding a subgraph with "<<subgraph->size()<<" nodes\n";
#endif

		if( Configure::REAL_POSITION_FILE != "" ){
			// label the real mapping for all contigs in this subgraph
			LabelContigsUsingMapping( subgraph );
		}

		m_isSpecialGraph = false;
		int numOfBC = 0;
		//if( subgraph->size() == 1330 )
		//	m_isSpecialGraph = true;

#ifdef DEBUG
		cerr<<"start a new subgraph\n";
		cerr<<"total contigs: "<<subgraph->size()<<endl;
#endif

		//printGraph = true;
		printGraph = false;
		int i = 0;
		for( list<Contig*>::iterator iter = subgraph->begin(); iter != subgraph->end(); iter++ ){
			i++;
			if( printGraph ){
			cerr<<i<<" contig\n";
			cerr<<(*iter)->GetName()<<"\tlength: "<<(*iter)->GetLength()<<" ";
			cerr<<"scaffold: \n"<<(*iter)->GetScaffoldString()<<endl;
			if( (*iter)->IsRepeat() )
				cerr<<"repeat\n";
			else
				cerr<<endl;
			}

			(*iter)->ClearUnassignedNodeIter( PLUS );
			(*iter)->ClearUnassignedNodeIter( MINUS );
			(*iter)->SetStartPosition( 1, 0 );
			(*iter)->SetIfCalculated( false );

			if( (*iter)->GetLength() > Configure::UPPERBOUND ) //Configure::LIB_MEAN + Configure::STD_TIMES * Configure::LIB_STD )
				numOfBC++;
			
			/*if( (*iter)->GetName() == "1056786" ){
				m_isSpecialGraph = false;
				}*/

#ifdef LOG
			// print contig
			//		fprintf( logFile, "%s\t%f\n", (*iter)->GetName().c_str(), (*iter)->GetLength() );
#endif
			
#ifdef LOG
//#ifdef DEBUG
			// print contig
			//fprintf( logFile, "%s\t%d\n", (*iter)->GetName().c_str(), (*iter)->GetLength() );

			// print 
			fprintf( logFile, "%s\tid is %d length is %f\n", (*iter)->GetName().c_str(), (*iter)->GetID(), (*iter)->GetLength() );
			list<PET*>::iterator edgeIter = (*iter)->GetLeftEdges()->begin();
			//cerr<<"left edges:\n";
			while( edgeIter != (*iter)->GetLeftEdges()->end() ){
				if( (*edgeIter)->IfInSubgraph() ){
					fprintf( logFile, "%s\n", (*edgeIter)->GetOriString().c_str() );
					//cerr<<"edge: "<<(*edgeIter)->GetOriString()<<endl;
				}
				else{
					//cerr<<"edge not in this subgraph: "<<(*edgeIter)->GetOriString()<<endl;
				}
				edgeIter++;
			}
			//cerr<<"right edges\n";
			edgeIter = (*iter)->GetRightEdges()->begin();
			while( edgeIter != (*iter)->GetRightEdges()->end() ){
				if( (*edgeIter)->IfInSubgraph() ){
					fprintf( logFile, "%s\n", (*edgeIter)->GetOriString().c_str() );
					//cerr<<"edge: "<<(*edgeIter)->GetOriString()<<endl;
				}
				else{
					//cerr<<"edge not in this subgraph: "<<(*edgeIter)->GetOriString()<<endl;
				}
				edgeIter++;
			}
			

			//fprintf( logFile, "%d\tlength: %.0f\n", (*iter)->GetID(), (*iter)->GetLength() );

#endif
		}

#ifdef DEBUG
		//cerr<<"finish graph"<<endl;
#endif

		
		if( m_isSpecialGraph ){
			GenerateDotFile( subgraph, 99 );

			/*for( list<Contig*>::iterator iter = subgraph->begin(); iter != subgraph->end(); iter++ )
			  cerr<<(*iter)->GetName()<<endl;*/

			//cerr<<"Total border contigs: "<<numOfBC<<endl;
		}
		

		// select starting contig
		/*Contig *startContig = NULL;
		int startContigOri;
		FindStartingContig( subgraph, numOfContig, numOfBorderContig, startContig, startContigOri );*/

		// scaffolding
		bool ifSucceed = false;
		ifIncreaseThreshold = false;		// record if this subgraph need to increase threshold
		
		int maxOfUnhappyEdges = -1;
		// create a scaffold containing startContig
		/*Clear();
		PartialScaffold *currentScaffold  = CreatScaffold( startContig, startContigOri );
		PartialScaffold *head = currentScaffold;
		GenUnassigndNodes( currentScaffold, true );*/

		// create initial scaffold
		//cerr<<"Creating initial scaffold\n";
		Clear();
		PartialScaffold *currentScaffold = CreateInitialScaffold( subgraph, numOfContig );
		//PartialScaffold *head = currentScaffold;
		//cerr<<"Finish creating initial scaffold\n";

		int numberOfScaffold = 0;
		struct timeval t_start,t_end;		// record the time for one run
		do{
			gettimeofday( &t_start, NULL );
			// starting from 0
			maxOfUnhappyEdges++;

			//cerr<<"Current maximum number of discordant edges is "<<maxOfUnhappyEdges<<endl;
#ifdef LOG
			fprintf( logFile, "Current maximum number of discordant edges is %d.\n", maxOfUnhappyEdges );
			fflush( logFile );
			//cout<<"Current maximum number of discordant edges is "<<maxOfUnhappyEdges<<endl;
#endif

			//pair<string, int> newPair( currentScaffold->GetScaffoldString(), 
			//	currentScaffold->GetNumOfUnhappyEdges() );
			//m_visitedTree->IfExist( currentScaffold->GetARString(),
			//		currentScaffold->GetDEString(), currentScaffold->GetNumOfUnhappyEdges() );
			//m_visitedScaffolds->insert( newPair );
			// check if current scaffold has more than required unhappy edges
			if( currentScaffold->GetNumOfUnhappyEdges() > maxOfUnhappyEdges ){
				maxOfUnhappyEdges++;
				continue;
			}

			numberOfScaffold = Scaffolding( currentScaffold, maxOfUnhappyEdges, subgraph->size() );
			fprintf( logFile, "Num of scaffolds is: %d\n",  numberOfScaffold );
#ifdef LOG
			//fprintf( logFile, "Total visited partial scaffolds are: %d\n",  m_visitedTree->GetTotalNum() );
#endif

			if( m_isSpecialGraph ){
				GenerateDotFile( subgraph, numOfSpecial++ );
			}

			//cerr<<"Total scaffolds: "<<numberOfScaffold<<endl;
			gettimeofday( &t_end, NULL );

			Clear();		// clear the map tree
			if( numberOfScaffold >= 1 ){
				// succeed
				ifSucceed = true;
			}
			else if( numberOfScaffold == Configure::TOO_MANY_PARTIAL_SCAFFOLDS ){
				int difTime = (t_end.tv_sec - t_start.tv_sec) + (t_end.tv_usec - t_start.tv_usec)/1000000.0;
				fprintf( logFile, "Total scaffolding time: %d seconds\n", 
					 difTime );
				fflush( logFile );
				m_totalNumberOfIncreasedSubgraph++;

				

#ifdef LOG
				fprintf( logFile, "Increase cluster threshold for current subgraph to %d.\n", minClusterSize + Configure::INCREASED_DELTA );
				fprintf( logFile, "Total nodes in subgraph: %d\n", (int) subgraph->size() );
				fflush( logFile );

				// output the graph information
				//GenerateDotFile( subgraph, m_graphIDIncreaseThreshold++ );
#endif
				string discardedEdges = "";
				// remove all clusters below new threshold
				int newThreshold = minClusterSize + Configure::INCREASED_DELTA;
				for( list<Contig*>::iterator iter = subgraph->begin(); iter != subgraph->end(); iter++ ){
					(*iter)->SetNotInSubgraph();
					(*iter)->SetExtensionOri( NEITHER );
					(*iter)->SetIfInAR( false );
					m_totalIncreasedNodes.insert( (*iter)->GetName() );

					list<PET*>::iterator edgeIter = (*iter)->GetLeftEdges()->begin();
					while( edgeIter != (*iter)->GetLeftEdges()->end() ){
						if( (*edgeIter)->IfInSubgraph() && (*edgeIter)->GetSize() < newThreshold ){
							// delete this edge
							if( !(*edgeIter)->IfVisited() )
								(*edgeIter)->VisitEdge();
							else{
								discardedEdges += (*edgeIter)->GetOriString() + "\n";
								m_graph->RemoveOneEdge();
							}

							edgeIter = (*iter)->GetLeftEdges()->erase( edgeIter );
						}
						else{
							(*edgeIter)->Initialize();
							edgeIter++;
						}
					}

					edgeIter = (*iter)->GetRightEdges()->begin();
					while( edgeIter != (*iter)->GetRightEdges()->end() ){
						if( (*edgeIter)->IfInSubgraph() && (*edgeIter)->GetSize() < newThreshold ){
							// delete this edge
							if( !(*edgeIter)->IfVisited() )
								(*edgeIter)->VisitEdge();
							else{
								m_graph->RemoveOneEdge();
								discardedEdges += (*edgeIter)->GetOriString() + "\n";
							}

							edgeIter = (*iter)->GetRightEdges()->erase( edgeIter );
						}
						else{
							(*edgeIter)->Initialize();
							edgeIter++;
						}
					}

					// check if this contig becomes a singleton
					if( !(*iter)->HasEdge() && !(*iter)->IsRepeat() )
						m_graph->SetAsSingleton( *iter );
	
				} // end traversing all contigs in subgraph

				// print all discarded edges
				fprintf( m_discardEdgeFile, "%s", discardedEdges.c_str() );

				// go back and regenerate subgraph
				ifSucceed = true;
				ifIncreaseThreshold = true;
				subgraph->clear();
				delete subgraph;
			}
			else{
				// did not find any scaffolds after searching the whole space
				//Clear();
				// set all contigs unassigned nodes iterator to null
				for( list<Contig*>::iterator iter = subgraph->begin(); iter != subgraph->end(); iter++ ){
					(*iter)->ClearUnassignedNodeIter( PLUS );
					(*iter)->ClearUnassignedNodeIter( MINUS );
					(*iter)->SetStartPosition( 1, 0 );
					(*iter)->SetIfCalculated( false );

					// initialize all edges
					list<PET*>::iterator edgeIter = (*iter)->GetLeftEdges()->begin();
					while( edgeIter != (*iter)->GetLeftEdges()->end() ){
						if( (*edgeIter)->IfInSubgraph() ){
							(*edgeIter)->SetDE( false );
							(*edgeIter)->SetHappy();
						}
						edgeIter++;
					}
					edgeIter = (*iter)->GetRightEdges()->begin();
					while( edgeIter != (*iter)->GetRightEdges()->end() ){
						if( (*edgeIter)->IfInSubgraph() ){
							(*edgeIter)->SetDE( false );
							(*edgeIter)->SetHappy();
						}
						edgeIter++;
					}
				}

				DeleteAllScaffolds( currentScaffold, NULL, true );
			        currentScaffold = CreateInitialScaffold( subgraph, numOfContig );
			}
			gettimeofday( &t_end, NULL );
			if( !ifSucceed && Configure::ABORT ){
				if( (t_end.tv_sec - t_start.tv_sec) + (t_end.tv_usec - t_start.tv_usec)/1000000.0 > 1800 ){
					// if running time is longer than 30 mins, then remind user if he wants to quit
					cout<<"The time required might be quite long, we suggest you to increase the value of \"cluster_threshold\" in configuration file"<<endl;
					cout<<"Quiting the program...Please run again with new parameters"<<endl;
					//delete currentScaffold;
					DeleteAllScaffolds( currentScaffold, NULL, true );
					return -1;
				}
			}
		} while( !ifSucceed );

		if( ifIncreaseThreshold )
		{ // go back to regenerate subgraph
			continue;
		}

		// update graph
#ifdef LOG
		/*
		if( numberOfScaffold > 1 )
			fprintf( logFile, "Generating %d scaffolds...\n", numberOfScaffold );
		else
			fprintf( logFile, "Generating %d scaffold...\n", numberOfScaffold );

		fflush( logFile );
		*/
#endif

		//cout<<numberOfScaffold<<endl;
		int numOfScaffold = 0;
		vector<vector<Contig*>*> *scafs = new vector<vector<Contig*>*>;
#ifdef LOG
		fprintf( logFile, "Starting separating scaffolds...\n" );
		fflush( logFile );
#endif
		//cerr<<"start separating scaffolds"<<endl;
		list<Contig*> *allRepeatCopies = new list<Contig*>;
		SeparateScaffolds( subgraph, numOfScaffold, scafs, currentScaffold, allRepeatCopies );
		//cerr<<"finish separating scaffolds"<<endl;
		// check connectivity
		//ScaffoldResult **result = new ScaffoldResult*[ numberOfScaffold ];   
		//GenerateResults( currentScaffold, result, numberOfScaffold );
#ifdef LOG
		fprintf( logFile, "Finish separating scaffolds\n" );
		fflush( logFile );
#endif

		// without checking connectivity
		//cerr<<"start generating results\n";
		//fprintf( logFile, "initializing results\n" );
		//fflush( logFile );
		ScaffoldResult **result = new ScaffoldResult*[ numOfScaffold ];   
		//fprintf( logFile, "finish initializing\n" );
		//fflush( logFile );
		GenerateResults( scafs, result );
		//cerr<<"finish generating results\n";

		for( int i = 0; i < (int) scafs->size(); i++ ){
			scafs->at( i )->clear();
			delete scafs->at( i );
		}
		scafs->clear();
		delete scafs;
	       

#ifdef LOG
		fprintf( logFile, "Updating related edges of assembled scaffolds...\n" );
		fflush( logFile );
#endif
		//cerr<<"updating graph\n";
		m_graph->UpdateGraph( subgraph, result, m_discardEdgeFile );
		//cerr<<"finish updating graph\n";

#ifdef LOG
		fprintf( logFile, "Finish current subgraph.\n\n" );
		fflush( logFile );
#endif

#ifdef DEBUG
		//PrintScaffoldResult( result, numberOfScaffold );
		//cout<<endl;
#endif

#ifdef DEBUG
		finishedScaffold++;
		cout<<"finish "<<itos( finishedScaffold )<<" subgraph"<<endl;
#endif
		subgraph->clear();
		//for( int i = 0; i < numberOfScaffold; i++ )
		for( int i = 0; i < numOfScaffold; i++ )
			delete result[ i ];
		delete[] result;
		// delete scaffold
		DeleteAllScaffolds( currentScaffold, NULL, false );    // problem is here

		// delete all repeat copies
		//cerr<<"start deleting..."<<endl;
		list<Contig*>::iterator rIter = allRepeatCopies->begin();
		while( rIter != allRepeatCopies->end() ){
			Contig *tempRepeat = *rIter;
			tempRepeat->DeleteEdgesForRepeat( LEFT );
			tempRepeat->DeleteEdgesForRepeat( RIGHT );
			delete tempRepeat;
			rIter = allRepeatCopies->erase( rIter );
			
		}
		//cerr<<"finish deleting..."<<endl;
		delete allRepeatCopies;
		
		
		delete subgraph;
	}
	
	//cout<<"Total abandend graph is "<<abandendGraph<<endl;
	//cout<<"Total abandend node is "<<abandendNode<<endl;
	return 1;
}

// generate the graph file for certain subgraph
void opera::GenerateDotFile( list<Contig*> *subgraph, int id ){
	string fileName = Configure::OUTPUT_FOLDER + "subgraph" + itos( id ) + ".dot";

	// write the dot file
	FILE *graphFile = fopen( fileName.c_str(), "w" );
	fprintf( graphFile, "graph gap%d {\nrankdir=LR;\noverlap=vpsc\n", id );

	//cerr<<"start printing contigs"<<endl;
	// print contigs
	/*for( list<Contig*>::iterator contigIter = subgraph->begin(); contigIter != subgraph->end(); contigIter++ ){
		Contig *currentContig = (*contigIter);
		fprintf( graphFile, "\"%s length: %.0f cov: %.1f\";\n", currentContig->GetName().c_str(), currentContig->GetLength(), currentContig->GetCov() );
		}*/
	//cerr<<"finish printing contigs\n";

	// print edges
	//cerr<<"start printing edges\n";
	for( list<Contig*>::iterator contigIter = subgraph->begin(); contigIter != subgraph->end(); contigIter++ ){
		Contig *currentContig = (*contigIter);

		PrintEdges( graphFile, currentContig, currentContig->GetLeftEdges(), LEFT, subgraph );
		PrintEdges( graphFile, currentContig, currentContig->GetRightEdges(), RIGHT, subgraph );
	}
	//cerr<<"finish printing edges\n";

	fprintf( graphFile, "}\n" );
	fclose( graphFile );

}

// print edges to file
// only print edges whose other contig ID is bigger than id of c
int opera::PrintEdges( FILE *file, Contig *c, list<PET*> *edges, int type, list<Contig*> *subgraph )
{
	string arrowtail = "";
	if( type == LEFT )
		arrowtail = "vee";
	else
		arrowtail = "crow";

	//cerr<<"Start handling edges: "<<edges->size()<<endl;
	int numOfEdge = 0;
	for( list<PET*>::iterator iter = edges->begin(); iter != edges->end(); iter++ ){
		Contig *otherContig = (*iter)->GetOtherContig( c );
		//cerr<<"other contig: "<<otherContig->GetName()<<"\t"<<"Length: "<<otherContig->GetLength()<<endl;

		if( c->GetID() > otherContig->GetID() )
			continue;

		int otherContigOri = (*iter)->GetOrientationOfContig( otherContig );
		int otherContigPos = (*iter)->GetPositionOfContig( otherContig );

		// do not output the contigs which is not in this subgraph
		if( !otherContig->isInSubgraph() )
			continue;

		string arrowhead = "";
		if( type == LEFT ){
			if( (otherContigPos == START && otherContigOri == PLUS ) 
			    || ( otherContigPos == END && otherContigOri == MINUS ) )
				arrowhead = "crow";
			else
				arrowhead = "vee";
		}
		else{
			if( ( otherContigPos == END && otherContigOri == PLUS ) 
			    || ( otherContigPos == START && otherContigOri == MINUS ) )
				arrowhead = "vee";
			else
				arrowhead = "crow";
		}

		//cerr<<"printing edges...\n";
		//cerr<<"first contig name: "<<c->GetName()<<" length: "<<c->GetLength()<<" cov: "<<c->GetCov()<<endl;
		//cerr<<"Other contig name: "<<otherContig->GetName()<<" length: "<<otherContig->GetLength()<<" cov: "<<otherContig->GetCov()<<endl;
		fprintf( file, "\"%s length: %.0f cov: %.1f\" -- \"%s length: %.0f cov: %.1f\" ", c->GetName().c_str(), c->GetLength(), c->GetCov(), 
			otherContig->GetName().c_str(), otherContig->GetLength(), otherContig->GetCov() );

		fprintf( file, 
			 "[fontsize=30, arrowsize=1, len=3.1, penwidth=1.80, arrowtail=%s, arrowhead=%s, style=dashed, dir=both, color=purple, label = \"dis: %d size: %d\"];\n", 
			 arrowtail.c_str(), arrowhead.c_str(), (*iter)->GetDis(), (*iter)->GetSize() );

		numOfEdge++;
	}
	//cerr<<"finish handling edges\n";

	return numOfEdge;
}

// generate unassigned contigs for initial scaffold
void opera::GenUnassigndNodesForInitialScaffold( StartPoint **allStart, int numOfStarts )
{
	// check if there is border contigs
	bool havingBorderContig = false;
	for( int i = numOfStarts - 1; i >= 0; i-- ){
		if( allStart[ i ]->GetContig()->GetLength() > Configure::UPPERBOUND ){ //Configure::LIB_MEAN + Configure::STD_TIMES * Configure::LIB_STD ){
			havingBorderContig = true;
		}
		allStart[ i ]->GetContig()->SetExtensionOri( NEITHER );
		allStart[ i ]->GetContig()->InitializeDistance();
		allStart[ i ]->GetContig()->ClearUnassignedNodeIter( allStart[ i ]->GetOri() );
		allStart[ i ]->GetContig()->SetStep( 1 );
	}
	
	
	if( havingBorderContig ){
	// add unassigned nodes
		//fprintf( logFile, "have border contigs\n" );
		int maxValue;
		if( Configure::QUICK_MODE ){
			if( Configure::TOP_VALUE < numOfStarts )
				maxValue = Configure::TOP_VALUE - 1;
			else
				maxValue = numOfStarts - 1 ;
		}
		else
			maxValue = numOfStarts - 1;

		for( int i = maxValue; i >= 0; i-- ){
			pair<Contig*, int> pair( allStart[ i ]->GetContig(), allStart[ i ]->GetOri() );
			m_unassignedNodes->push_front( pair );
			pair.first->SetUnassignedNodeIter( pair.second, m_unassignedNodes->begin() );
			//InsertContigIntoUnassignedNodeSet( pair );
		}
	}
	else{
		// sort according to the distance

		for( list< pair<Contig*, int> >::iterator iter = m_unassignedNodes->begin(); iter != m_unassignedNodes->end(); iter++ ){
			(*iter).first->SetExtensionOri( NEITHER );
			(*iter).first->InitializeDistance();
		}
	
		// check new unassigned nodes
		list< pair<Contig*, int> > *possibleContig = new list< pair<Contig*, int> >;	// int is the extension direction
		set<Contig*> *newUnassignedNodes = new set<Contig*>;
		set<PET*> *visitedEdge = new set<PET*>;
		
		// deal with all contigs in active region
		Contig *firstContig = allStart[ 0 ]->GetContig();
		for( int i = 0; i < numOfStarts; i++ ){
			if( allStart[ i ]->GetContig()->IsRepeat()){
				firstContig = allStart[ i ]->GetContig();
				break;
			}
		}
		
		//fprintf( logFile, "the start contigs is: %s\n", firstContig->GetName().c_str() );
		firstContig->SetStep( 1 );
		firstContig->SetDistance( 0, 1 );
		firstContig->SetOriInExtension( PLUS );
		
		firstContig->SetExtensionOri( BOTH );
		int ext = BOTH;
		pair<Contig*, int> p( firstContig, ext );
		possibleContig->push_back( p );
		newUnassignedNodes->insert( firstContig );
		
		// start traversing
		bool isStart = true;
		bool ifRecalculate = true;
		while( !possibleContig->empty() ){
			// get the first contig
			pair<Contig*, int> pair = *possibleContig->begin();
			Contig *currentContig = pair.first;
			int ext = pair.second;
			possibleContig->pop_front();
			//fprintf( logFile, "traverse from %s\n", currentContig->GetName().c_str() );

			if( currentContig->IsRepeat() ){
				// do not extend from repeat
				continue;
			}
			
			// traverse edges
			if( ext == RIGHT )
				TraverseEdges( currentContig, currentContig->GetRightEdges(), possibleContig,
					       newUnassignedNodes, visitedEdge, RIGHT, isStart, ifRecalculate );
			else if( ext == LEFT )
				TraverseEdges( currentContig, currentContig->GetLeftEdges(), possibleContig,
					       newUnassignedNodes, visitedEdge, LEFT, isStart, ifRecalculate );
			else if( ext == BOTH ){
				TraverseEdges( currentContig, currentContig->GetLeftEdges(), possibleContig,
					       newUnassignedNodes, visitedEdge, LEFT, isStart, ifRecalculate );
				TraverseEdges( currentContig, currentContig->GetRightEdges(), possibleContig,
					       newUnassignedNodes, visitedEdge, RIGHT, isStart, ifRecalculate );
			}
		}
		
		// create the array
		StartPoint **allStarts;
		int num = 0;
		allStarts = new StartPoint*[ newUnassignedNodes->size() * 2 ];
		for( set<Contig*>::iterator iter = newUnassignedNodes->begin(); iter != newUnassignedNodes->end(); iter++ ){
			allStarts[ num++ ] = new StartPoint( *iter, PLUS );
			allStarts[ num++ ] = new StartPoint( *iter, MINUS );
		}

		//fprintf( logFile, "total unassigned nodes: %d\n", num );
		
		sort( allStarts, allStarts + num, SortStartPointsUsingDis );

		// add unassigned nodes
		int maxValue;
		if( Configure::QUICK_MODE ){
			if( Configure::TOP_VALUE < num )
				maxValue = Configure::TOP_VALUE - 1;
			else
				maxValue = num - 1 ;
		}
		else
			maxValue = num - 1;
		for( int i = maxValue; i >= 0; i-- ){
			pair<Contig*, int> pair( allStarts[ i ]->GetContig(), allStarts[ i ]->GetOri() );
			m_unassignedNodes->push_front( pair );
			pair.first->SetUnassignedNodeIter( pair.second, m_unassignedNodes->begin() );
			//fprintf( logFile, "%s pos: %f\n", allStarts[ i ]->GetContig()->GetName().c_str(), allStarts[ i ]->GetContig()->GetMidPointDistance() );
			//InsertContigIntoUnassignedNodeSet( pair );
		}

		//?// Change added
		for( int i = 0; i < num; i++ ){
			delete allStarts[i];
		}
		delete[] allStarts;



		newUnassignedNodes->clear();   delete newUnassignedNodes;
		possibleContig->clear();  delete possibleContig;
		visitedEdge->clear();  delete visitedEdge;
	}
	   
}

// create the initial scaffold containing all nodes as unassigned nodes
PartialScaffold* opera::CreateInitialScaffold( list<Contig*> *subgraph, int numOfContig )
{
	StartPoint **allStart;
	int num = 0;
	//int type;

	allStart = new StartPoint*[ numOfContig * 2 ];
	//type = ALL;

	list<Contig*>::const_iterator iter;
	for( iter = subgraph->begin(); iter != subgraph->end(); iter++ ){
		// do not use repeat as the start contig in scaffolds
		if( (*iter)->IsRepeat() )
			continue;

		allStart[ num++ ] = new StartPoint( *iter, PLUS );
		allStart[ num++ ] = new StartPoint( *iter, MINUS );
	}

	sort( allStart, allStart + num, SortStartPoints );

	/*if( m_isSpecialGraph ){
		cerr<<"After sorting:\n";
		for( int i = 0; i < num; i++ ){
			cerr<<allStart[ i ]->GetContig()->GetName()<<"\t"<<allStart[ i ]->GetContig()->GetLength()<<"\t"<<allStart[ i ]->GetNumOfUnhappyEdges()<<"\t";
			if( allStart[ i ]->GetOri() == PLUS )
				cerr<<" + \n";
			else
				cerr<<" - \n";
			//cerr<<allStart[ i ]->GetContig()->GetScaffoldString()<<endl<<endl;
			
		}
	}
	*/

	PartialScaffold *result = new PartialScaffold;

	// generate unassigned node
	//GenUnassigndNodesForInitialScaffold( allStart, numOfContig*2 );
	GenUnassigndNodesForInitialScaffold( allStart, num );

	//cerr<<"Total unassigned nodes: "<<m_unassignedNodes->size()<<endl;

	for( int i = 0; i < num; i++ )
		delete allStart[ i ];

	delete[] allStart;

		
	result->SetFirstScaffold();

	return result;
}

// separate the scaffolds into independent ones
void opera::SeparateScaffolds( list<Contig*> *subgraph, int &numberOfScaf, vector<vector<Contig*>*> *scafs, PartialScaffold *currentScaffold,
			       list<Contig*> *allRepeatCopies )
{
#ifdef DEBUG
	fprintf( logFile, "Total contig in subgraph: %d\n", subgraph->size() );
	fflush( logFile );
#endif
	//cout<<"Total contig in subgraph: "<<subgraph->size()<<endl;
	
	// collect all the nodes in scaffolds (the repeats which are actually used)
	list<Contig*> newSubgraph;
	// traverse scaffolds
	PartialScaffold *tempPS = currentScaffold;

#ifdef DEBUG
	fprintf( logFile, "Adding contigs into graph\n" );
	fflush( logFile );
#endif
	//cerr<<"scaffold:\n";
	//cerr<<"contigs in subgraph:\n";
	while( tempPS->GetParent() != NULL ){
		if( tempPS->IfEmptyScaffold() ){
			// it is a new scaffold
			
			// start a new scaffold
			if( tempPS->IfMannualEmptyScaffold() )
				tempPS = tempPS->GetParent();
		}

		Contig *addedContig = tempPS->GetAddedContigInAR();
		newSubgraph.push_back( addedContig );
		if( addedContig->IsRepeat() ){
			allRepeatCopies->push_back( addedContig );
		}
		//cerr<<addedContig->GetName()<<endl;
		
		/*
		cerr<<addedContig->GetName()<<"\t";
		if( addedContig->GetOri() == PLUS )
			cerr<<"+\t";
		else
			cerr<<"-\t";
		*/

		tempPS = tempPS->GetParent();
	}
	//cerr<<endl<<endl;
	
#ifdef DEBUG
	fprintf( logFile, "Finish adding contigs into graph\n" );
	fflush( logFile );
#endif

	// traverse graph and label the scaffold ID for each contig
	set<Contig*> *visitedContig = new set<Contig*>;
	list<Contig*> *queue = new list<Contig*>;
	int scafID = -1;
	
	// save unhappy edges information
	set<PET*> *unhappyEdges = new set<PET*>;

	/*for( list<Contig*>::iterator graphIter = newSubgraph.begin(); graphIter != newSubgraph.end(); graphIter++ ){
		Contig *newContig = *graphIter;
		if( !newContig->IsRepeat() && visitedContig->find( newContig ) == visitedContig->end() ){
			// this node is unique and has not been visited, start traverse the graph to find the connected component
			queue->push_back( newContig );
			visitedContig->insert( newContig );
			scafID++;

			while( !queue->empty() ){
				Contig *currentContig = queue->front();
				queue->pop_front();
				
				currentContig->SetScafID( scafID );
				currentContig->SetScaffoldID( scafID );

				// do not extend from repeat
				if( currentContig->IsRepeat() )
					continue;

				// get all connected contigs
				list<PET*> *edges = new list<PET*>;
				edges->insert( edges->end(), currentContig->GetLeftEdges()->begin(), currentContig->GetLeftEdges()->end() );
				edges->insert( edges->end(), currentContig->GetRightEdges()->begin(), currentContig->GetRightEdges()->end() );
				for( list<PET*>::iterator edgeIter=edges->begin(); edgeIter != edges->end(); edgeIter++ ){
					if( !(*edgeIter)->isInSubgraph() )
						continue;

					if( (*edgeIter)->IfUnhappy() ){
						if( !Configure::SIMULTANEOUSLY_HANDLE_LIBRARIES )
							unhappyEdges->insert( *edgeIter );
						else{
							for( list<PET*>::iterator singlePetIter = (*edgeIter)->GetClustersOfSuperCluster()->begin();
							     singlePetIter != (*edgeIter)->GetClustersOfSuperCluster()->end(); singlePetIter++ )
								unhappyEdges->insert( *singlePetIter );
						}
						continue;
					}

					Contig *nextContig = (*edgeIter)->GetOtherContig( currentContig );
					if( nextContig->IsRepeat() ){
						// nextContig is repeat, get the real occurrence
						cerr<<"unique contig is :"<<currentContig->GetName()<<endl;
						cerr<<(*edgeIter)->GetOriString()<<endl;
						cerr<<"start contig: "<<(*edgeIter)->GetStartContig()->GetName()<<"\tend contig: "<<(*edgeIter)->GetEndContig()->GetName()<<endl;
						if( (*edgeIter)->GetRepetitiveEdge() != NULL )
							nextContig = (*edgeIter)->GetRepetitiveEdge()->GetOtherContig( currentContig );
						else
							nextContig = (*edgeIter)->GetUniqueEdge()->GetOtherContig( currentContig );
						//nextContig = (*edgeIter)->GetOtherContig( currentContig );
					}

					if( visitedContig->find( nextContig ) == visitedContig->end() ){
						// if nextContig has not been visited before, add this contig into this scaffold
						queue->push_back( nextContig );
						visitedContig->insert( nextContig );
					}
				}
				edges->clear();
				delete edges;
			}
		}
	}
	*/

#ifdef DEBUG
	fprintf( logFile, "Traversing graph\n" );
	fflush( logFile );
#endif
	//cerr<<"Traversing graph\n";
	// new traversing to handle uniqe edges (single direction edge)
	vector<set<Contig*>* > scafSets;
	for( list<Contig*>::iterator graphIter = newSubgraph.begin(); graphIter != newSubgraph.end(); graphIter++ ){
		Contig *newContig = *graphIter;
		//cerr<<"node: "<<newContig->GetName()<<endl;
		if( !newContig->IsRepeat() && visitedContig->find( newContig ) == visitedContig->end() ){
			// this node is unique and has not been visited, start traverse the graph to find the connected component
			//cerr<<"one graph\n";
			set<Contig*> *newSet = new set<Contig*>;
			scafSets.push_back( newSet );

			queue->push_back( newContig );
			//visitedContig->insert( newContig );
			//scafID = scafSets.size();
			scafID = scafSets.size() - 1;

			while( !queue->empty() ){
				Contig *currentContig = queue->front();
				queue->pop_front();

#ifdef DEBUG
	fprintf( logFile, "visiting %s\n", currentContig->GetName().c_str() );
	fflush( logFile );
#endif
				
				// if this contig is visited before, skip it
				if( visitedContig->find( currentContig ) != visitedContig->end() ){
					continue;
				}
			       
				currentContig->SetScafID( scafID );
				currentContig->SetScaffoldID( scafID );
				currentContig->SetConnectedSetPos( scafID );
				visitedContig->insert( currentContig );
				newSet->insert( currentContig );
				//cerr<<"visit: "<<currentContig->GetName()<<"\t"<<"connected set pos is "<<currentContig->GetConnectedSetPos()<<endl;
				//cerr<<"number of nodes in set: "<<newSet->size()<<endl;

				// do not extend from repeat
				if( currentContig->IsRepeat() )
					continue;

				// get all connected contigs
				list<PET*> *edges = new list<PET*>;
				edges->insert( edges->end(), currentContig->GetLeftEdges()->begin(), currentContig->GetLeftEdges()->end() );
				edges->insert( edges->end(), currentContig->GetRightEdges()->begin(), currentContig->GetRightEdges()->end() );
				for( list<PET*>::iterator edgeIter=edges->begin(); edgeIter != edges->end(); edgeIter++ ){
					if( !(*edgeIter)->isInSubgraph() )
						continue;

					if( (*edgeIter)->IfUnhappy() ){
						if( !Configure::SIMULTANEOUSLY_HANDLE_LIBRARIES )
							unhappyEdges->insert( *edgeIter );
						else{
							for( list<PET*>::iterator singlePetIter = (*edgeIter)->GetClustersOfSuperCluster()->begin();
							     singlePetIter != (*edgeIter)->GetClustersOfSuperCluster()->end(); singlePetIter++ )
								unhappyEdges->insert( *singlePetIter );
						}
						continue;
					}

					Contig *nextContig = (*edgeIter)->GetOtherContig( currentContig );
					if( nextContig->IsRepeat() ){
						// nextContig is repeat, get the real occurrence
						//cerr<<"unique contig is :"<<currentContig->GetName()<<endl;
						//cerr<<(*edgeIter)->GetOriString()<<endl;
						//cerr<<"start contig: "<<(*edgeIter)->GetStartContig()->GetName()<<"\tend contig: "<<(*edgeIter)->GetEndContig()->GetName()<<endl;
						if( (*edgeIter)->GetRepetitiveEdge() != NULL ){
							//cerr<<"get repetitive edge\n";
							nextContig = (*edgeIter)->GetRepetitiveEdge()->GetOtherContig( currentContig );
						}
						else{
							//cerr<<"get unique edge\n";
							nextContig = (*edgeIter)->GetUniqueEdge()->GetOtherContig( currentContig );
						}

						//cerr<<"finish hanlding repeat"<<endl;
						//nextContig = (*edgeIter)->GetOtherContig( currentContig );
					}

					if( visitedContig->find( nextContig ) == visitedContig->end() ){
						// if nextContig has not been visited before, add this contig into this scaffold
						//cerr<<"push back to queue\n";
						queue->push_back( nextContig );
						//visitedContig->insert( nextContig );
					}
					else{
						//cerr<<"next contig: "<<nextContig->GetName()<<endl;
						if( nextContig->GetConnectedSetPos() != scafID ){
							// combine two sets
							//cerr<<"combine two sets\n";
							//cerr<<"old set id is "<<nextContig->GetConnectedSetPos()<<endl;
							set<Contig*> *combinedSet = scafSets.at( nextContig->GetConnectedSetPos() );
							//cerr<<"number of nodes in combined set: "<<combinedSet->size()<<endl;
							for( set<Contig*>::iterator tempSetIter = newSet->begin(); tempSetIter != newSet->end(); tempSetIter++ ){
								//cerr<<"add: "<<(*tempSetIter)->GetName()<<endl;
								(*tempSetIter)->SetScafID( nextContig->GetConnectedSetPos() );
								(*tempSetIter)->SetScaffoldID( nextContig->GetConnectedSetPos() );
								(*tempSetIter)->SetConnectedSetPos( nextContig->GetConnectedSetPos() );
								combinedSet->insert( (*tempSetIter) );
							}
							newSet->clear();
							newSet = combinedSet;
							scafID = nextContig->GetConnectedSetPos();
							//cerr<<"number of nodes in combined set: "<<combinedSet->size()<<endl;
						}
					}
				}
				edges->clear();
				delete edges;
			}
		}
	}

	//cerr<<"finish Traversing\n";
#ifdef DEBUG
	fprintf( logFile, "Finish traversing graph\n" );
	fflush( logFile );
#endif

	scafID = -1;
	for( int i = 0; i < (int) scafSets.size(); i++ ){
		if( scafSets.at( i )->size() > 0 ){
			scafID++;
			for( set<Contig*>::iterator tempIter = scafSets.at( i )->begin(); tempIter != scafSets.at( i )->end(); tempIter++ ){
				(*tempIter)->SetScafID( scafID );
				(*tempIter)->SetScaffoldID( scafID );
				(*tempIter)->SetConnectedSetPos( scafID );
				//cerr<<(*tempIter)->GetName()<<"\t"<<(*tempIter)->GetScafID()<<endl;
			}
		}
		scafSets.at( i )->clear();
		delete scafSets.at( i );
	}
	scafSets.clear();

#ifdef DEBUG
	fprintf( logFile, "Finish reindexing\n" );
	fflush( logFile );
#endif
	//cerr<<"finish reindexing\n";

	for( set<PET*>::iterator iter = unhappyEdges->begin(); iter != unhappyEdges->end(); iter++ ){
		//cout<<(*iter)->GetOriString()<<endl;
		m_unhappyEdgeString->push_back( (*iter)->GetOriString() );
	}
	
	//fprintf( logFile, "finish separate scaffolds\n" );
	//fflush( logFile );
	//cerr<<"finish separating scaffolds\n";
	
	visitedContig->clear();
	delete visitedContig;
	queue->clear();
	delete queue;
	unhappyEdges->clear();
	delete unhappyEdges;
	
#ifdef DEBUG
	fprintf( logFile, "start saving scaffolds\n" );
	fflush( logFile );
#endif
	//cerr<<"start saving scaffolds\n";
	// traverse the scaffolds and save each scaffold
	numberOfScaf = scafID + 1;
	//cout<<numberOfScaf<<endl<<endl;
	// initialize the scaffolds
	//cerr<<"number of scaffold: "<<numberOfScaf<<endl;
	for( int i = 0; i < numberOfScaf; i++ )
		scafs->push_back( new vector<Contig*> );

	// generate scaffolds
	tempPS = currentScaffold;
	int contigOrderID = 0;
	//int numberOfGaps = 0;

	int numOfContig = 0;
	//Contig *lastAddedContig = NULL;
	while( tempPS->GetParent() != NULL ){
		if( tempPS->IfEmptyScaffold() ){
			// it is a new scaffold
			//cerr<<"start a new scaffold\n";
			
			// start a new scaffold
			if( tempPS->IfMannualEmptyScaffold() )
				tempPS = tempPS->GetParent();

			contigOrderID = 0;
		}

		if( tempPS->GetPosOfUnassignedNode() > 0 ){
			if( tempPS->GetPosOfUnassignedNode() <= 10 ){
				m_numOfTopUnassignedNodes++;
			}
			//cerr<<tempPS->GetPosOfUnassignedNode()<<"\t"<<tempPS->GetAddedContigInAR()->GetStep()<<"\t"<<tempPS->GetAddedContigInAR()->GetName()<<endl;
			m_totalPartialScaffold++;
		}

		Contig *addedContig = tempPS->GetAddedContigInAR();
			//addedContig->SetScaffoldPos( contigOrderID );
		contigOrderID++;
	
		// alternatives
		//cerr<<"set scaffold position for contig: "<<addedContig->GetName()<<endl;
		//cerr<<"scaffold id: "<<addedContig->GetScafID()<<endl;
		
		addedContig->SetScaffoldPos( scafs->at( addedContig->GetScafID() )->size() );
		scafs->at( addedContig->GetScafID() )->push_back( addedContig );
		//cerr<<"finish setting\n";

		tempPS = tempPS->GetParent();
		numOfContig++;
		//lastAddedContig = addedContig;
		//cout<<addedContig->GetName()<<"\t"<<addedContig->GetOri()<<endl;
	}
	
	//cout<<endl;
	//cout<<"Total contig in scaffold: "<<numOfContig<<endl<<endl;

#ifdef DEBUG
	// print all scaffold
	/*for( int i = 0; i < scafs->size(); i++ ){
		cerr<<"scaffold "<<i<<":"<<endl;
		for( int j = 0; j < scafs->at( i )->size(); j++ ){
			cerr<<scafs->at( i )->at( j )->GetName()<<"\t";
			if( scafs->at( i )->at( j )->GetOri() == PLUS )
				cerr<<"+\t";
			else
				cerr<<"-\t";
			
		}
		cerr<<endl<<endl;
	}
	*/

	/*for( int i = 0; i < scafs->size(); i++ ){
		cerr<<"scaffold "<<i<<":"<<endl;
		for( int j = 0; j < scafs->at( i )->size(); j++ ){
			if( scafs->at( i )->at( j )->GetName() == "1455" ){
				// print the edge
				for( list<PET*>::iterator tempEdgeIter = scafs->at( i )->at( j )->GetRightEdges()->begin(); 
				     tempEdgeIter != scafs->at( i )->at( j )->GetRightEdges()->end(); tempEdgeIter++ ){
					cerr<<(*tempEdgeIter)->GetOriString()<<endl;
					cerr<<(*tempEdgeIter)<<endl;
				}
			}
			
		}
		cerr<<endl<<endl;
	}
	*/
#endif
}

// delete all scaffold
inline void opera::DeleteAllScaffolds( PartialScaffold *p, PartialScaffold *head, bool clearRepeat ){
	while( p != head ){
		PartialScaffold *temp = p;
		
		if( clearRepeat ){
			// delete occurrences of all repeats
			Contig *currentContig = p->GetAddedContigInAR();
			//if( temp->IfAddedContigIsRepeat() && currentContig != NULL ){
			if( currentContig != NULL && currentContig->IsRepeat() ){
				currentContig->DeleteEdgesForRepeat( currentContig->GetLeftEdges() );
				currentContig->DeleteEdgesForRepeat( currentContig->GetRightEdges() );
				delete currentContig;
			}
		}	 

		p = p->GetParent();
		delete temp;
	}
}

bool SortStartPoints( const StartPoint *s1, const StartPoint *s2 ){
	if( s1->GetNumOfUnhappyEdges() == s2->GetNumOfUnhappyEdges() ){
		return s1->GetContig()->GetLength() > s2->GetContig()->GetLength();
		//return s1->GetContig()->GetDistance() > s2->GetContig()->GetDistance();
	}
	else
		return s1->GetNumOfUnhappyEdges() < s2->GetNumOfUnhappyEdges();
}

bool SortStartPointsUsingDis( const StartPoint *s1, const StartPoint *s2 ){
	return s1->GetContig()->GetMidPointDistance() < s2->GetContig()->GetMidPointDistance();
}

// select the starting contig
void opera::FindStartingContig( list<Contig*> *subgraph, int numOfContig, int numOfBorderContig,
							   Contig *&startContig, int &startContigOri ){
	StartPoint **allStart;
	int num = 0;
	int type;
	if( numOfBorderContig != 0 ){
		allStart = new StartPoint*[ numOfBorderContig * 2 ];
		type = BORDER;
	}
	else{
		allStart = new StartPoint*[ numOfContig * 2 ];
		type = ALL;
	}

	list<Contig*>::const_iterator iter;
	for( iter = subgraph->begin(); iter != subgraph->end(); iter++ ){
		if( ( (type == BORDER ) && ( (*iter)->IsBorderContig() ) ) 
			|| (type == ALL ) ){
				allStart[ num++ ] = new StartPoint( *iter, PLUS );
				allStart[ num++ ] = new StartPoint( *iter, MINUS );
		}
	}

	sort( allStart, allStart + num, SortStartPoints );
	startContig = allStart[ 0 ]->GetContig();
	startContigOri = allStart[ 0 ]->GetOri();
	for( int i = 0; i < num; i++ )
		delete allStart[ i ];

	delete[] allStart;
}	

// create a scaffold containing a start contig
PartialScaffold* opera::CreatScaffold( Contig *c, int ori ){

	PartialScaffold *result = new PartialScaffold;
		
	// add a contig to active region
	c->SetOri( ori );
	c->SetStartPosition( 1, 0 );
	result->AddNode( c );
	m_activeRegion->push_back( c );
	
	// add edges to dangling edges
	if( ori == PLUS ){
		// add all right edges
		AddHappyDEToFirstScaffold( c->GetRightEdges(), result );
		AddUnhappyDEToFirstScaffold( c->GetLeftEdges(), result, c );
	}
	else{
		// add all left edges
		AddHappyDEToFirstScaffold( c->GetLeftEdges(), result );
		AddUnhappyDEToFirstScaffold( c->GetRightEdges(), result, c );
	}

	GenARString( result );
	GenUDEString( result );

	// generate unassigned node
		
	return result;
}

// add happy dangling edges to first scaffolds
void opera::AddHappyDEToFirstScaffold( list<PET*> *edges, PartialScaffold *p ){
	list<PET*>::const_iterator iter = edges->begin();
	for( iter = edges->begin(); iter != edges->end(); iter++ ){
		if( (*iter)->isInSubgraph() ){
			p->AddAddedHappyDE( *iter );
			m_happyDanglingEdges->push_front( *iter );
			(*iter)->SetHappyDEIter( m_happyDanglingEdges->begin() );
		}
	}
}

// add unhappy dangling edges to first scaffolds
// c is the first contig
void opera::AddUnhappyDEToFirstScaffold( list<PET*> *edges, PartialScaffold *p, Contig *c ){
	list<PET*>::const_iterator iter = edges->begin();
	for( iter = edges->begin(); iter != edges->end(); iter++ ){
		if( (*iter)->isInSubgraph() ){
			p->AddAddedUnhappyDE( *iter, (*iter)->GetOtherContig( c ) );
			//cout<<"unhappy due to first contig:"<<(*iter)->GetOriString()<<endl;
			m_unhappyDanglingEdges->push_front( *iter );
			(*iter)->SetUnhappyDEIter( m_unhappyDanglingEdges->begin() );
		}
	}
}

bool UDESort( PET *&p1, PET *&p2 ){
	return p1->GetID() < p2->GetID();
}

// generate unhappy dangling edges string
void opera::GenUDEString( PartialScaffold *p ){
	if( m_unhappyDanglingEdges->empty() ){
		p->SetUnhappyDEString( "" );
		return;
	}

	// sort unhappy dangling edges
	m_unhappyDanglingEdges->sort( UDESort );

	string result = "";
	for( list<PET*>::const_iterator iter = m_unhappyDanglingEdges->begin();
		iter != m_unhappyDanglingEdges->end(); iter++ ){
			result += itos((*iter)->GetID()) + "\t";
	}

	p->SetUnhappyDEString( result );
}

// print AR string
string opera::PrintARStringOfName(){
	string result = "";

	for( list<Contig*>::const_iterator iter = m_activeRegion->begin(); 
		iter != m_activeRegion->end(); iter++ ){
			result += (*iter)->GetName() + "\t";
			if( (*iter)->GetOri() == PLUS )
				result += "+\t";
			else
				result += "-\t";
	}

	return result;
}

// generate active region string
void opera::GenARString( PartialScaffold *p ){
	string result = "";

	for( list<Contig*>::const_iterator iter = m_activeRegion->begin(); 
		iter != m_activeRegion->end(); iter++ ){
			result += itos( (*iter)->GetID() ) + "\t";
			if( (*iter)->GetOri() == PLUS )
				result += "+\t";
			else
				result += "-\t";
	}

	p->SetARString( result );
}

// scaffolding
// return -1 if failed
// otherwise, return the number of scaffold
int opera::Scaffolding( PartialScaffold *&currentScaffold, int maxOfUnhappyEdges, int numberOfContigs ){
	int numberOfScaffold = 1;

	// initialize the tree saving all visited scaffolds
	m_visitedTree->Clear();

	// while it is not root partial scaffold or it is not finished
	//while( currentScaffold->GetParent() != NULL || !currentScaffold->IfBreaked() ){
	int numOfTrial = 0;
	while( currentScaffold->GetParent() != NULL || !m_unassignedNodes->empty() ){	
		//numOfTrial++;
#ifdef DEBUG
		//cerr<<currentScaffold->GetAddedContigInAR()->GetName()<<endl;
		cerr<<"previous scaffolds: \n";
		// traverse scaffolds
		PartialScaffold *tempPS = currentScaffold;

		cerr<<"scaffold:\n";
		while( tempPS->GetParent() != NULL ){
			if( tempPS->IfEmptyScaffold() ){
				// start a new scaffold
				if( tempPS->IfMannualEmptyScaffold() )
					tempPS = tempPS->GetParent();
			}
			
			Contig *addedContig = tempPS->GetAddedContigInAR();
			
			cerr<<addedContig->GetName()<<"\t";
			if( addedContig->GetOri() == PLUS )
				cerr<<"+\t";
			else
				cerr<<"-\t";
			
			tempPS = tempPS->GetParent();
		}
		cerr<<endl<<endl;
	
		
		cerr<<"current active region:\n";
		for( list<Contig*>::iterator cIter = m_activeRegion->begin(); cIter != m_activeRegion->end(); cIter++ ){
			cerr<<(*cIter)->GetName()<<" ";
			if( (*cIter)->GetOri() == PLUS )
				cerr<<"+ ";
			else
				cerr<<"- ";
		}
		cerr<<"\n";
		
		/*cerr<<"happy dangling edges:\n";
		for( list<PET*>::iterator dIter = m_happyDanglingEdges->begin(); dIter != m_happyDanglingEdges->end(); dIter++ ){
			cerr<<(*dIter)->GetOriString()<<endl;
			if( !(*dIter)->IsDE() || (*dIter)->IfUnhappy() )
				cerr<<"not a happy dangling edge\n";
		}
		cerr<<endl;
		cerr<<"unhappy dangling edges:\n";
		for( list<PET*>::iterator dIter = m_unhappyDanglingEdges->begin(); dIter != m_unhappyDanglingEdges->end(); dIter++ ){
			cerr<<(*dIter)->GetOriString()<<endl;
		}
		cerr<<endl;
		*/
		
		
#endif
		//if( m_visitedTree->GetTotalNum() > Configure::MAX_NUMBER_OF_PARTIAL_SCAFFOLDS )
		if( numOfTrial > Configure::MAX_NUMBER_OF_PARTIAL_SCAFFOLDS )
		{
			// Trace back the partial scaffold
			while( currentScaffold->GetParent() != NULL ){ //!currentScaffold->IsFirstScaffold() ){
				TraceBack( currentScaffold );
			}
			
			//cerr<<"too many trials: "<<numOfTrial<<endl;
			fprintf( logFile, "total trial: %d\n", numOfTrial );

			return Configure::TOO_MANY_PARTIAL_SCAFFOLDS;
		}

		if( currentScaffold->IsFirstScaffold() || currentScaffold->IfEmptyScaffold() || !m_activeRegion->empty() ){
			// select an unassigned node and add to the end
			numOfTrial++;
			//cerr<<numOfTrial<<endl;

			Contig *contig = NULL;
			int contigOri = 0;
			
#ifdef DEBUG
			if( !m_unassignedNodes->empty() ){
				cerr<<"There are "<<m_unassignedNodes->size()<<" unassinged nodes\n";
				cerr<<"The first one is "<<m_unassignedNodes->front().first->GetName()<<endl;
				for( list< pair<Contig*, int> >::iterator tempUnassIter = m_unassignedNodes->begin(); tempUnassIter != m_unassignedNodes->end(); tempUnassIter++ ){
					cerr<<(*tempUnassIter).first->GetName()<<"\t";
					if( (*tempUnassIter).second == PLUS )
						cerr<<"+\t";
					else
						cerr<<"-\t";
				       
				}
				cerr<<endl<<endl;
			}
			else{
				cerr<<"there is no unassigned nodes\n";
			}
#endif
					

			// active region is not empty
			if( !m_unassignedNodes->empty() ){
				// if there still unused unassigned nodes
				// try to add a contig
				//cerr<<"Try to add node\n";
				if( !AddContigToScaffold( currentScaffold, maxOfUnhappyEdges, contig, contigOri ) ){
					// if adding is not successful, trace back
					// new scaffold is not valid, need to trace back
#ifdef DEBUG
					cerr<<"Adding failed\n";
#endif
					TraceBack( currentScaffold );
					if( currentScaffold->IfEmptyScaffold() ){
						//numberOfScaffold--;
					}
				}
				else{
					// adding successful, check if new scaffold has been visited yet or not
					// check if new scaffold exist and insert if not
					/*if( m_isSpecialGraph ){
						cerr<<"Adding succeed!\n";
						cerr<<currentScaffold->GetAddedContigInAR()->GetName()<<endl;
						}*/
#ifdef DEBUG
					cerr<<"Adding correct\n";
#endif
					bool ifExist = m_visitedTree->IfExist( currentScaffold->GetARString(),
						currentScaffold->GetDEString(), currentScaffold->GetNumOfUnhappyEdges() );
					if( ifExist ){
						// no need to continue, trace back
						TraceBack( currentScaffold );
						if( currentScaffold->IfEmptyScaffold() ){
#ifdef SPLIT
							cout<<"minus 1 for number of scaffolds"<<endl;
							if( currentScaffold->IfMannualEmptyScaffold() )
								cout<<"return from a manual splited one"<<endl;
							else
								cout<<"retrun from a normal splited one"<<endl;
#endif
							//numberOfScaffold--;
						}	
					}

#ifdef SPLIT
					string logString = "\nmax unhappy edges: " + itos( maxOfUnhappyEdges ) + "\n";
					logString += "new Scaffold\n" + *(currentScaffold->GetARString()) + "\n" + *(currentScaffold->GetDEString()) + "\n";
					fprintf( logFile, logString.c_str() );
					//cerr<<"finish adding\n";
#endif
				}
			}
			else{
				// if there is no unassigned nodes
				// trace back
				TraceBack( currentScaffold );
				continue;
			}	
		}
		else{
			//cerr<<"active region is empty\n";
			// active region is empty
			if( m_unhappyDanglingEdges->empty() ){
				// if there is no unhappy dangling edges, find a solution
				//cerr<<"no unhappy dangling edges, find a solution.\n";
#ifdef SPLIT
				fprintf( logFile, "return after adding contig\n" );
#endif
				//cerr<<"Find the results\n";
				return numberOfScaffold;
			}
			else{
				// if there is some unhappy dangling edges, start a new scaffold
				//cerr<<"Still has some unhappy dangling edges, start a new scaffold.\n";
				currentScaffold->SetIfEmptyScaffold( true );
				currentScaffold->Break();
				//numberOfScaffold++;
#ifdef SPLIT
				cout<<"split a new scaffold "<<numberOfScaffold<<endl;
				cout<<"last contig in last scaffold is "<<currentScaffold->GetAddedContigInAR()->GetName()<<endl;
#endif
				// need to start a new scaffold
				// find and sort all the unassigned nodes for this new scaffold
				CreateNewScaffold();

				/*StartPoint *sp = FindStartOfNewScaffold();
				sp->GetContig()->SetStartPosition( 1, 0 );
				// create a scaffold
				if( !AddContigToScaffold( currentScaffold, maxOfUnhappyEdges, 
							sp->GetContig(), sp->GetOri() ) ){
						// cannot start new scaffold, trace back
						TraceBack( currentScaffold );	// track back to end of last scaffold
						TraceBack( currentScaffold );	// change an end
						numberOfScaffold--;
#ifdef SPLIT
						cout<<"trace back from AddContigToScaffold"<<endl;
#endif
						if( currentScaffold->IfEmptyScaffold() )
							numberOfScaffold--;
				}
				else{
					// create new scaffold successfully
					if( m_visitedTree->IfExist( currentScaffold->GetARString(),
								currentScaffold->GetDEString(), currentScaffold->GetNumOfUnhappyEdges() ) ){
									// if this new scaffold already exist
							TraceBack( currentScaffold );	// track back to end of last scaffold
							TraceBack( currentScaffold );	// change an end
							numberOfScaffold--;
#ifdef SPLIT
						cout<<"trace back from being been visited before"<<endl;
#endif
						if( currentScaffold->IfEmptyScaffold() )
							numberOfScaffold--;
					}
				}
				delete sp;*/
			}
		}
	}
	//cerr<<"Total trial: "<<numOfTrial<<endl;
	fprintf( logFile, "total trial: %d\n", numOfTrial );

	return -1;
}

// create a new scaffold for all disconnected nodes
void opera::CreateNewScaffold(){
	// get all the disconnected nodes
	list<Contig*> *subgraph = new list<Contig*>;

	// traverse
	set<Contig*> *visitedContigs = new set<Contig*>;
	list<Contig*> *possibleContigs = new list<Contig*>;
		// start with unhappy dangling edges
	for( list<PET*>::iterator iter = m_unhappyDanglingEdges->begin();
		iter != m_unhappyDanglingEdges->end(); iter++ ){
			Contig *contig = (*iter)->GetUnusedContig();
			//cerr<<"unused contig: "<<contig->GetName()<<endl;
			if( visitedContigs->insert( contig ).second ){
				// it is the first time for this contig, put into the list
				possibleContigs->push_back( contig );
				subgraph->push_back( contig );
			}
	}

	while( !possibleContigs->empty() ){
		Contig *contig = *(possibleContigs->begin());
		possibleContigs->pop_front();

		CheckEdgesForUnusedContigs( subgraph, contig, contig->GetLeftEdges(), visitedContigs, possibleContigs );
		CheckEdgesForUnusedContigs( subgraph, contig, contig->GetRightEdges(), visitedContigs, possibleContigs );
	}

	// create scaffold
	PartialScaffold *newScaf;
	newScaf = CreateInitialScaffold( subgraph, subgraph->size() );
	delete newScaf;

	// delete variables
	visitedContigs->clear();
	delete visitedContigs;
	possibleContigs->clear();
	delete possibleContigs;
	subgraph->clear();
	delete subgraph;
}

// Break scaffold
/*void opera::BreakScaffold( PartialScaffold *&p, int maxUnhappyEdges ){
	PartialScaffold *emptySP = new PartialScaffold();
	emptySP->SetNumOfUnhappyEdges( p->GetNumOfUnhappyEdges() );
	emptySP->SetParent( p );
	emptySP->Break();
	emptySP->SetIfEmptyScaffold( true );
	emptySP->SetIfMannualEmptyScaffold( true );

	int increasedUnhappyEdges = m_happyDanglingEdges->size();
	// set all happy dangling edges as unhappy dangling edges
	if( increasedUnhappyEdges + p->GetNumOfUnhappyEdges() > maxUnhappyEdges ){
		delete emptySP;
		return;
	}

	list<PET*>::iterator edgeIter = m_happyDanglingEdges->begin();
	while( edgeIter != m_happyDanglingEdges->end() ){
		PET *edge = *edgeIter;
		edgeIter = m_happyDanglingEdges->erase( edgeIter );
		m_unhappyDanglingEdges->push_front( edge );
		edge->SetUnhappyDEIter( m_unhappyDanglingEdges->begin() );
		emptySP->AddRemovedHappyDE( edge );
		emptySP->AddAddedUnhappyDE( edge, edge->GetOtherContig( edge->FindARContig() ) );
	}

	list<Contig*>::iterator contigIter = m_activeRegion->begin();
	while( contigIter != m_activeRegion->end() ){
		Contig *contig = *contigIter;
		contigIter = m_activeRegion->erase( contigIter );
		contig->SetIfInAR( false );
		emptySP->AddRemovedContigInAR( contig );
		
		// check if currentContig connects to active region
		//if( !CheckActiveRegionConnectivity( contig ) ){
		//	TraceBack( emptySP );
		//	return;
		//	}
	}

	StartPoint *sp = FindStartOfNewScaffold();
	// create a scaffold
	if( !AddContigToScaffold( emptySP, maxUnhappyEdges, sp->GetContig(), sp->GetOri() ) ){
		// cannot start new scaffold
		//PartialScaffold *parent = emptySP->GetParent(); 
		TraceBack( emptySP );
		TraceBack( emptySP );
		delete sp;
		return;
	}

	delete sp;

	p = emptySP;
}
*/


// find start contig of new scaffold
// return start point
StartPoint* opera::FindStartOfNewScaffold(){
	list<StartPoint*> *borderContigs = new list<StartPoint*>;	// border contigs
	list<StartPoint*> *allContigs = new list<StartPoint*>;		// non border contigs
	StartPoint *result = NULL;

	// traverse
	set<Contig*> *visitedContigs = new set<Contig*>;
	list<Contig*> *possibleContigs = new list<Contig*>;
		// start with unhappy dangling edges
	for( list<PET*>::iterator iter = m_unhappyDanglingEdges->begin();
		iter != m_unhappyDanglingEdges->end(); iter++ ){
			Contig *contig = (*iter)->GetUnusedContig();
			if( visitedContigs->insert( contig ).second ){
				// it is the first time for this contig, put into the list
				AddToList( borderContigs, allContigs, contig, PLUS );
				AddToList( borderContigs, allContigs, contig, MINUS );
				possibleContigs->push_back( contig );
			}
	}

	while( !possibleContigs->empty() ){
		Contig *contig = *(possibleContigs->begin());
		possibleContigs->pop_front();

		CheckEdgesForStartPoint( borderContigs, allContigs, contig, contig->GetLeftEdges(), 
			visitedContigs, possibleContigs );
		CheckEdgesForStartPoint( borderContigs, allContigs, contig, contig->GetRightEdges(), 
			visitedContigs, possibleContigs );
	}

	if( !borderContigs->empty() ){
		borderContigs->sort( SortStartPoints );
		result = *(borderContigs->begin());
		borderContigs->pop_front();
	}
	else{
		allContigs->sort( SortStartPoints );
		result = *(allContigs->begin());
		allContigs->pop_front();
	}

	// delete variables
	list<StartPoint*>::iterator iter = borderContigs->begin();
	while( iter != borderContigs->end() ){
		StartPoint *sp = *iter;
		iter = borderContigs->erase( iter );
		delete sp;
	}
	delete borderContigs;
	iter = allContigs->begin();
	while( iter != allContigs->end() ){
		StartPoint *sp = *iter;
		iter = allContigs->erase( iter );
		delete sp;
	}
	delete allContigs;
	delete visitedContigs;
	delete possibleContigs;

	return result;
}

// check edges to find all unused contigs
void opera::CheckEdgesForUnusedContigs( list<Contig*> *subgraph, Contig *c, list<PET*> *edges, 
					set<Contig*> *visitedContigs, list<Contig*> *possibleContigs ){
	for( list<PET*>::iterator iter = edges->begin(); iter != edges->end(); iter++ ){
		PET *edge = *iter;
		if( edge->isInSubgraph() && !edge->IsDE() && !edge->IfUnhappy() ){
			Contig *otherContig = edge->GetOtherContig( c );
			if( otherContig->isInSubgraph() && !otherContig->IsRepeat() && visitedContigs->insert( otherContig ).second ){
				// it is the first time for this contig, put into the list
				possibleContigs->push_back( otherContig );
				subgraph->push_back( otherContig );
			}
		}
	}
}

// check edges to find start point
void opera::CheckEdgesForStartPoint( list<StartPoint*> *borderContigs, list<StartPoint*> *allContigs, 
									Contig *c, list<PET*> *edges, set<Contig*> *visitedContigs,
									list<Contig*> *possibleContigs ){
	for( list<PET*>::iterator iter = edges->begin(); iter != edges->end(); iter++ ){
		PET *edge = *iter;
		if( edge->isInSubgraph() && !edge->IsDE() && !edge->IfUnhappy() ){
			Contig *otherContig = edge->GetOtherContig( c );
			if( otherContig->isInSubgraph() && visitedContigs->insert( otherContig ).second ){
				// it is the first time for this contig, put into the list
				AddToList( borderContigs, allContigs, otherContig, PLUS );
				AddToList( borderContigs, allContigs, otherContig, MINUS );
				possibleContigs->push_back( otherContig );
			}
		}
	}
}

// add connected contigs to corresponding list
void opera::AddToList( list<StartPoint*> *borderContigs, list<StartPoint*> *allContigs, Contig *c, int ori ){
	StartPoint *newSP = new StartPoint( c, ori );
	newSP->SetNumOfUnhappyEdges( c->GetNumOfValidLeftEdges( ori ) );
	if( c->IsBorderContig() )
		borderContigs->push_back( newSP );
	else
		allContigs->push_back( newSP );
}

// generate unassigned contig set
// isStart means if this is the start of the scaffold
void opera::GenUnassigndNodes( PartialScaffold *p, bool isStart, bool ifRecalculate ){
	// clear unassigned nodes possible orientation
	for( list< pair<Contig*, int> >::iterator iter = m_unassignedNodes->begin(); iter != m_unassignedNodes->end(); iter++ ){
		(*iter).first->SetExtensionOri( NEITHER );
		(*iter).first->InitializeDistance();
	}

	// check new unassigned nodes
	list< pair<Contig*, int> > *possibleContig = new list< pair<Contig*, int> >;	// int is the extension direction
	set<Contig*> *newUnassignedNodes = new set<Contig*>;
	set<PET*> *visitedEdge = new set<PET*>;

	// deal with all contigs in active region
	//if( m_isSpecialGraph)
	//	cerr<<"Active region edges:\n";
	for( list<Contig*>::iterator iter = m_activeRegion->begin(); iter != m_activeRegion->end(); iter++ ){
		Contig *contig = *iter;
		if( contig->IsRepeat() ){
			// do not use repeat to extend the scaffold
			continue;
		}

		list<PET*> *edges;
		if( (*iter)->GetOri() == PLUS ){
			edges = (*iter)->GetRightEdges();
			if( isStart || ifRecalculate )
				contig->SetOriInExtension( PLUS );
		}
		else{
			edges = (*iter)->GetLeftEdges();
			if( isStart || ifRecalculate )
				contig->SetOriInExtension( MINUS );
		}

		for( list<PET*>::iterator edgeIter = edges->begin(); edgeIter != edges->end(); edgeIter++ ){
			if( !(*edgeIter)->IfUnhappy() && (*edgeIter)->isInSubgraph() && (*edgeIter)->IsDE() ){
				// only consider happy edge
				PET *edge = *edgeIter;
				visitedEdge->insert( edge );
				int ori = edge->GetOrientationOfContig( contig );
				Contig *otherContig = edge->GetOtherContig( contig );
				int otherOri = edge->GetOrientationOfContig( otherContig );
				otherContig->SetStep( 1 );
				
				// calculate the distance of this contig to the end of active region
				int disToTheEndOfActiveRegion = 0;
				for( list<Contig*>::reverse_iterator riter = m_activeRegion->rbegin(); riter != m_activeRegion->rend(); riter++ ){
					if( *riter != contig )
						//disToTheEndOfActiveRegion += (*riter)->GetLength() + (*riter)->GetDistance();
						disToTheEndOfActiveRegion += (*riter)->GetLength() + (*riter)->GetDistance();
					else
						break;
				}

				// set the distance if this is for the first scaffold
				if( isStart || ifRecalculate ){
					//int tempDis = 0;
					if( otherContig->IfCalculated() ){
						// already have a distance, check if need to combine the distances
						//otherContig->AddDistance( edge->GetDis() );
						otherContig->AddDistance( edge->GetDis() - disToTheEndOfActiveRegion, 1 );
					}
					else{
						//otherContig->SetDistance( edge->GetDis() );
						otherContig->SetDistance( edge->GetDis() - disToTheEndOfActiveRegion, 1 );
						//otherContig->SetStep( 0 );
						//fprintf( logFile, "%s  dis: %f\n", otherContig->GetName().c_str(), otherContig->GetDistance() );
						//fflush( logFile );
						if( ori == contig->GetOriInExtension() )
							otherContig->SetOriInExtension( otherOri );
						else
							otherContig->SetOriInExtension( GetOppositeOri( otherOri ) );
					}

					if( otherContig->GetName() == "629503"  ){
						//cerr<<"Traverse edge from "<<(*iter)->GetName()<<"\t";
						//cerr<<"to "<<otherContig->GetName()<<" dis: "<<otherContig->GetDistance()<<"\tvisited time: "<<otherContig->GetVisitedTime()<<endl;
						fprintf( logFile, "Traverse edge from %s\t to %s dis: %f\tvisited time: %d\n",
							 (*iter)->GetName().c_str(), otherContig->GetName().c_str(), otherContig->GetDistance(), 
							 otherContig->GetVisitedTime() );
					}
				}
				int ext;
				if( !otherContig->IsBorderContig() ){
					// it is not border contig, extend both direction
					otherContig->SetExtensionOri( BOTH );
					ext = BOTH;
				}
				else{
					// it is border contig, only extend to one direction
					if( ( ori == contig->GetOri() && otherOri == PLUS ) ||
						( ori != contig->GetOri() && otherOri == MINUS ) ){
							otherContig->SetExtensionOri( LEFT );
							ext = LEFT;
					}
					else{
						otherContig->SetExtensionOri( RIGHT );
						ext = RIGHT;
					}
				}
				pair<Contig*, int> p( otherContig, ext );
				possibleContig->push_back( p );
				newUnassignedNodes->insert( otherContig );
			}
		}
	}

	//if( m_isSpecialGraph)
	//	cerr<<"Traversing edges\n";

	// start traversing
	while( !possibleContig->empty() ){
		// get the first contig
		pair<Contig*, int> pair = *possibleContig->begin();
		Contig *currentContig = pair.first;
		int ext = pair.second;
		possibleContig->pop_front();

		if( currentContig->IsRepeat() ){
			// do not extend from repeat
			continue;
		}

		if( Configure::QUICK_MODE ){
			// only consider contigs within several steps
			if( currentContig->GetStep() >= Configure::STEP )
				continue;
		}

		// traverse edges
		if( ext == RIGHT )
			TraverseEdges( currentContig, currentContig->GetRightEdges(), possibleContig,
				       newUnassignedNodes, visitedEdge, RIGHT, isStart, ifRecalculate );
		else if( ext == LEFT )
			TraverseEdges( currentContig, currentContig->GetLeftEdges(), possibleContig,
				       newUnassignedNodes, visitedEdge, LEFT, isStart, ifRecalculate );
		else if( ext == BOTH ){
			TraverseEdges( currentContig, currentContig->GetLeftEdges(), possibleContig,
				       newUnassignedNodes, visitedEdge, LEFT, isStart, ifRecalculate );
			TraverseEdges( currentContig, currentContig->GetRightEdges(), possibleContig,
				       newUnassignedNodes, visitedEdge, RIGHT, isStart, ifRecalculate );
		}
	}

	// find removed unassigned nodes and update corresponding variables
	list< pair<Contig*, int> >::iterator pairIter = m_unassignedNodes->begin();
	while( pairIter != m_unassignedNodes->end() ){
		if( (*pairIter).first->GetExtensionOri() == NEITHER ){
			// it is the removed unassigned nodes
			p->AddRemovedUnassignedNodes( (*pairIter).first, (*pairIter).second );
			(*pairIter).first->ClearUnassignedNodeIter( (*pairIter).second );
			pairIter = m_unassignedNodes->erase( (*pairIter).first->GetUnassignedNodeIter( (*pairIter).second ) );
		}
		else
			pairIter++;
	}
		
	// find added unassigned nodes
	for( set<Contig*>::iterator iter = newUnassignedNodes->begin(); iter != newUnassignedNodes->end(); iter++ ){
		if( !(*iter)->IsUnassignedNode( PLUS ) ){
			// not an unassigned nodes before
			p->AddAddedUnassignedNodes( *iter, PLUS );
			pair<Contig*, int> pair1( *iter, PLUS );
			InsertContigIntoUnassignedNodeSet( pair1 );
		}
		if( !(*iter)->IsUnassignedNode( MINUS ) ){
			p->AddAddedUnassignedNodes( *iter, MINUS );
			pair<Contig*, int> pair2( *iter, MINUS );
			InsertContigIntoUnassignedNodeSet( pair2 );
		}
	}

	//cerr<<"recalculating...\n";
	// if need to re-sort all unassigned nodes
	if( ifRecalculate){
		// clear the previous unassigned nodes set
		list< pair<Contig*, int> >::iterator unassignedIter = m_unassignedNodes->begin();
		while( unassignedIter != m_unassignedNodes->end() ){
			(*unassignedIter).first->ClearUnassignedNodeIter( PLUS );
			(*unassignedIter).first->ClearUnassignedNodeIter( MINUS );
			unassignedIter++;
		}

		m_unassignedNodes->clear();

		// add all new unassigned nodes into the set
		// find added unassigned nodes
		bool ifPrint = false;
		int num = 0;
		for( set<Contig*>::iterator iter = newUnassignedNodes->begin(); iter != newUnassignedNodes->end(); iter++ ){
			if( Configure::QUICK_MODE ){
				num++;
				if( num > Configure::TOP_VALUE )
					break;
			}

			// not an unassigned nodes before
			pair<Contig*, int> pair1( *iter, PLUS );
			InsertContigIntoUnassignedNodeSet( pair1 );
			
			pair<Contig*, int> pair2( *iter, MINUS );
			InsertContigIntoUnassignedNodeSet( pair2 );
		}
		if( ifPrint )
			cerr<<endl<<endl;
	}

	//cerr<<"First unassigned nodes is : "<<(m_unassignedNodes->begin())->first->GetName()<<endl;

	newUnassignedNodes->clear();   delete newUnassignedNodes;
	possibleContig->clear();  delete possibleContig;
	visitedEdge->clear();  delete visitedEdge;
}

// insert a contig into unassigned nodes set according to distance
void opera::InsertContigIntoUnassignedNodeSet( pair<Contig*, int> pa ){
	list< pair<Contig*, int> >::iterator iter;

	// if the distance of the new contig is less than -2000, than insert it into the end of the list
	//if( pa.first->GetDistance() < -2000 ){
	/*if( pa.first->GetMidPointDistance() < -2000 ){ // use 6 standard deviation
		iter = m_unassignedNodes->end();
		m_unassignedNodes->insert( iter, pa );
		iter--;
		pa.first->SetUnassignedNodeIter( pa.second, iter );
		return;
	}
	*/

	// for all the contigs with less than 5 steps, compare them together; 
	// never compare the contigs with less than 5 steps with more than 5 steps;
	for( iter = m_unassignedNodes->begin(); iter != m_unassignedNodes->end(); iter++ ){
		//bool lessThan10 = (*iter).first->GetStep() > 10 && pa.first->GetStep() <= 10;
		//bool inSameRegion = ((*iter).first->GetStep() <= 10 && pa.first->GetStep() <= 10 ) 
		//	|| ((*iter).first->GetStep() > 10 && pa.first->GetStep() > 10);
		//bool lessDis = (*iter).first->GetDistance() >= pa.first->GetDistance();
		bool lessDis = abs( (*iter).first->GetMidPointDistance() ) >= abs( pa.first->GetMidPointDistance() );
		//if( lessThan10 || (inSameRegion && lessDis) ){
		if( lessDis ){
			// insert the new contig into this position
			m_unassignedNodes->insert( iter, pa );
			iter--;
			pa.first->SetUnassignedNodeIter( pa.second, iter );
			return;
		}
	}

	// insert into the end of the list
	m_unassignedNodes->insert( iter, pa );
	iter--;
	pa.first->SetUnassignedNodeIter( pa.second, iter );
	return;
}

// traverse edge and find unassigned node
void opera::TraverseEdges( Contig *c, list<PET*> *edges, list< pair<Contig*, int> > *&possibleContig,
			   set<Contig*> *&unassignedNodes, set<PET*> *&visitedEdges, int direction, bool isStart, bool ifRecalculate ){
	list<PET*>::iterator iter;
	for( iter = edges->begin(); iter != edges->end(); iter++ ){
		if( !(*iter)->isInSubgraph() || (*iter)->IfUnhappy() || (*iter)->IsDE() ){
			continue;
		}

		/*if( visitedEdges->find(*iter) != visitedEdges->end() )
			continue;
		*/

		else{
			visitedEdges->insert(*iter);
			Contig *otherContig = (*iter)->GetOtherContig( c );
			int ori = (*iter)->GetOrientationOfContig( c );
			int otherPos = (*iter)->GetPositionOfContig( otherContig );
			int otherOri = (*iter)->GetOrientationOfContig( otherContig );
			int ext;

			if( isStart || ifRecalculate ){
				// if the distance of the other contig has not been calcuated yet
				// calculate the distance
				double dis = 0;
				if( c->GetOriInExtension() == PLUS ){
					if( direction == RIGHT ){
						// plus distance
						dis = c->GetDistance() + c->GetLength() + (*iter)->GetDis();
					}
					else{
						// minus distance
						dis = c->GetDistance() - (*iter)->GetDis() - otherContig->GetLength();
					}
				}
				else{
					if( direction == RIGHT ){
						// minus distance
						dis = c->GetDistance() - (*iter)->GetDis() - otherContig->GetLength();
					}
					else{
						// plus distance
						dis = c->GetDistance() + c->GetLength() + (*iter)->GetDis();
					}
				}

				if( otherContig->IfCalculated() ){
					// add this distance
					otherContig->AddDistance( dis, c->GetStep() + 2 );
					//if( otherContig->GetName() == "629503" ){
					if( m_isSpecialGraph ){
						//cerr<<"Reach not the first time:\n";
						//cerr<<"dis: "<<dis<<endl;
						//cerr<<"Traverse edge from "<<c->GetName()<<"\t"<<"dis: "<<c->GetDistance()<<"\t";
						//cerr<<"to "<<otherContig->GetName()<<" dis: "<<otherContig->GetDistance()<<"\tvisited time: "<<otherContig->GetVisitedTime()<<endl;
						fprintf( logFile, "Reach not the first time:\n" );
						fprintf( logFile, "dis: %f\n", dis );
						fprintf( logFile, "Traverse edge from %s\t to %s dis: %f\tvisited time: %d\n",
							 c->GetName().c_str(), otherContig->GetName().c_str(), otherContig->GetDistance(), 
							 otherContig->GetVisitedTime() );
					}
				}
				else{
					// set the orientation of other contig
					if( ori == c->GetOriInExtension() ){
						otherContig->SetOriInExtension( otherOri );
					}
					else{
						otherContig->SetOriInExtension( GetOppositeOri( otherOri ) );
					}

					otherContig->SetStep( c->GetStep() + 1 );
					otherContig->SetDistance( dis, c->GetStep() + 2 );
					//if( otherContig->GetName() == "629503" ){
					if( m_isSpecialGraph ){
						//cerr<<"Reach the first time:\n";
						//cerr<<"Traverse edge from "<<c->GetName()<<"\t"<<"dis: "<<c->GetDistance()<<"\t";
						//cerr<<"to "<<otherContig->GetName()<<" dis: "<<otherContig->GetDistance()<<"\tvisited time: "<<otherContig->GetVisitedTime()<<endl;
						fprintf( logFile, "Reach the first time:\n" );
						fprintf( logFile, "Traverse edge from %s\t to %s dis: %f\tvisited time: %d\n",
							 c->GetName().c_str(), otherContig->GetName().c_str(), otherContig->GetDistance(), 
							 otherContig->GetVisitedTime() );
					}
				}


				//fprintf( logFile, "%s  dis: %f\n", otherContig->GetName().c_str(), otherContig->GetDistance() );
				//fflush( logFile );
			}

			int previousExt = otherContig->GetExtensionOri();
			if( !otherContig->IsBorderContig() ){
				// short contig, just add and extend to both orientation
				otherContig->SetExtensionOri( BOTH );
				ext = BOTH;
			}
			else{
				// border contig, need consider orientation
				if( ( otherPos == START && otherOri == MINUS ) 
				|| ( otherPos == END && otherOri == PLUS ) ){
					// new contig should extend to left
					otherContig->SetExtensionOri( LEFT );
					ext = LEFT;
				}
				else{
					otherContig->SetExtensionOri( RIGHT );
					ext = RIGHT;
				}
			}
			int newExt = otherContig->GetExtensionOri();
			
			// if this contig needs to be extended to a new direction, put it in the possibleContig list
			if( previousExt != newExt ){
				otherContig->SetStep( c->GetStep() + 1 );
				pair<Contig*, int> p(otherContig, ext);
				possibleContig->push_back( p );
				unassignedNodes->insert( otherContig );
			}
		}
	}
}


// add one contig to partial scaffold
// if the new scaffold is not valid, return false
// if newContig is NULL, select one contig from unassigned node list
//   or just add newContig to the end of scaffold
bool opera::AddContigToScaffold( PartialScaffold *&p, int maxUnhappyEdge, Contig *newContig, 
								int newContigOri ){

	if( Configure::REAL_POSITION_FILE != "" ){
		pair<Contig*, int> firstPair = *m_unassignedNodes->begin();
		Contig *addContig = firstPair.first;

		// check if there is any confliction between the real position and the trial
		if( p->GetParent() == NULL ){
			// it is the root scaffold, check if the trial is the first contig in reference
			if( addContig->GetSubgraphReferenceIndex() != 0 && addContig->GetSubgraphReferenceIndex() != m_subgraphSize - 1){
				fprintf( logFile, "WARNING: The start contig %s is not the first contig in this subgraph\n", addContig->GetName().c_str() );
				fprintf( logFile, "index: %d, length: %f, total contigs: %d\n", addContig->GetSubgraphReferenceIndex(), 
					 addContig->GetLength(), m_subgraphSize );

				//m_isSpecialGraph = true;
			}
		}
		else{
			// check if the adjacent contigs are really adjacent in reference
			if( !m_activeRegion->empty() ){
				Contig *endContig = m_activeRegion->back();	
				fprintf( logFile, "Active region: %s\n", PrintARStringOfName().c_str() );
				if( abs( endContig->GetSubgraphReferenceIndex() - addContig->GetSubgraphReferenceIndex() ) != 1 ){
					fprintf( logFile, "WARNING: The next contig %s is not the adjacent one of the previous contig %s in this graph\n", 
						 addContig->GetName().c_str(), endContig->GetName().c_str() );
					fprintf( logFile, "Distance of this contig is %d.\n", (int) addContig->GetDistance() );
					//cerr<<"WARNING: The next contig "<<addContig->GetName()<<" is not the adjacent one of the previous contig"
					//    <<endContig->GetName()<<" in this graph\n";
					//cerr<<"Distance of this contig is "<<addContig->GetDistance()<<endl;

					//m_isSpecialGraph = true;
				}
			}
			
		}
	}

	// generate a new partial scaffold
	PartialScaffold *newScaffold = new PartialScaffold( p );
	newScaffold->SetParent( p );
	newScaffold->SetNumOfUnhappyEdges( p->GetNumOfUnhappyEdges() );
	p = newScaffold;

	// get the added contig
	//cerr<<"Add contig\n";
	bool needToReturn = false;
	if( newContig == NULL ){
		p->GetParent()->IncreasePosOfUnassignedNode();
		pair<Contig*, int> firstPair = *m_unassignedNodes->begin();

		newContig = firstPair.first;
		newContigOri = firstPair.second;
		m_unassignedNodes->erase( newContig->GetUnassignedNodeIter( newContigOri ) );
		newContig->ClearUnassignedNodeIter( newContigOri );

		//cerr<<"left edges: "<<newContig->GetLeftEdges()->size()<<endl;
		//cerr<<"right edges: "<<newContig->GetRightEdges()->size()<<endl;

#ifdef DEBUG 
		cerr<<"Adding:\t";
		if( newContig->IsRepeat() )
			cerr<<"repeat\t";

		cerr<<newContig->GetName()<<" ";
		if( newContigOri == PLUS )
			cerr<<"+\n";
		else
			cerr<<"-\n";
		
#endif

		// record the added contig is not repeat
		p->SetTypeOfAddedContig( false );

		if( newContig->IsRepeat() ){
			//cerr<<"is repeat\n";
			// create another contig to represent this repeat
			Contig *repeatContig = newContig;

			// record the added contig is repeat
			p->SetTypeOfAddedContig( true );

			// if the newly added contig is the same as the last contig in active region, do not add this contig
			Contig *lastAddedContig = p->GetParent()->GetAddedContigInAR();
			if( lastAddedContig != NULL ){
				if( lastAddedContig->IsRepeat() ){
					if( lastAddedContig->GetOriginalRepeat() == repeatContig 
					    && lastAddedContig->GetOri() == newContigOri ){
						needToReturn = true;
						//return false;
					}
				}
				else if( lastAddedContig->GetLastSubcontigs()->m_contigName == repeatContig->GetName() 
					 && lastAddedContig->GetLastSubcontigs()->m_ori == newContigOri ){
					needToReturn = true;
					//return false;
				}
			}

			newContig = new Contig( repeatContig->GetName(), repeatContig->GetLength(), repeatContig->GetCov() );
			newContig->SetID( repeatContig->GetID() );
			newContig->SetIfRepeat( true );
			newContig->SetOriginalRepeat( repeatContig );
			newContig->ClearUnassignedNodeIter( newContigOri );

			int otherOri = GetOppositeOri( newContigOri );
			if( repeatContig->IsUnassignedNode( otherOri ) ){
				newContig->SetUnassignedNodeIter( otherOri, repeatContig->GetUnassignedNodeIter( otherOri ) );
			}
			else{
				newContig->ClearUnassignedNodeIter( otherOri );
			}
			
			//cerr<<"create a new contig for this repeat\n";
								  
		}

		newContig->SetOri( newContigOri );

		// if this contig is a repeat, create subcontigs for it
		if( newContig->IsRepeat() ){
			newContig->CreateSubContigForRepeat();
		}
		
		/*fprintf( logFile, "try to add %s\t", newContig->GetName().c_str() );
		if( newContigOri == PLUS )
			fprintf( logFile, "+\n" );
		else
			fprintf( logFile, "-\n" );
		*/
		

		if( !m_activeRegion->empty() ){
			Contig *endContig = m_activeRegion->back();		// the last contig in active region
			newContig->SetStartPosition( endContig->GetStartPosition() + endContig->GetLength() - endContig->GetRemovedDis(), 0 );
		}
		else{
			newContig->SetStartPosition( 1, 0 );
		}
		//m_unassignedNodes->erase( newContig->GetUnassignedNodeIter( newContigOri ) );
		//newContig->ClearUnassignedNodeIter( newContigOri );

		p->SavePreviousUnassignedNodesOrder( m_unassignedNodes );

		/*if( m_isSpecialGraph ){
			cerr<<"Add the contig: "<<newContig->GetName()<<"\t"<<newContig->GetLength()<<"\t";
			if( newContigOri == PLUS )
				cerr<<" + "<<endl;
			else
				cerr<<" - "<<endl;
				}*/
		

		//cout<<(*p->GetARString())<<endl;
		//cout<<"adding: "<<newContig->GetName()<<"\t"<<newContigOri<<endl;
		//cout<<newContig->GetScaffoldString()<<endl<<endl;
	}
	else{
		newContig->SetOri( newContigOri );
		newContig->SetStartPosition( 1, 0 );
	}

	if( needToReturn ){
		// fix me: delete new contig if it is a repeat
		if( newContig->IsRepeat() ){
			delete newContig;
		}

		return false;
	}

	//cerr<<"check all edges of new contig\n";
	// check all edges of new Contig
	// FIXME: if an edge with repeat is satisfied, create a unique edge
	// fixed
	if( !newContig->IsRepeat() ){
		CheckEdgesOfNewContig( newScaffold, newContig, newContig->GetLeftEdges(), newContigOri );
		//cerr<<"discordant edges: "<<newScaffold->GetNumOfUnhappyEdges()<<endl;
		CheckEdgesOfNewContig( newScaffold, newContig, newContig->GetRightEdges(), newContigOri );
		//cerr<<"discordant edges: "<<newScaffold->GetNumOfUnhappyEdges()<<endl;
	}
	else{
		CheckEdgesOfNewContig( newScaffold, newContig, newContig->GetOriginalRepeat()->GetLeftEdges(), newContigOri );
		//cerr<<"discordant edges: "<<newScaffold->GetNumOfUnhappyEdges()<<endl;
		CheckEdgesOfNewContig( newScaffold, newContig, newContig->GetOriginalRepeat()->GetRightEdges(), newContigOri );
		//cerr<<"discordant edges: "<<newScaffold->GetNumOfUnhappyEdges()<<endl;
	}
	
	//cerr<<"add v to partial scaffold\n";
	// add v to partial scaffold
	newScaffold->AddNode( newContig );
	m_activeRegion->push_back( newContig );
	if( newScaffold->GetParent() != NULL && !newScaffold->GetParent()->IfEmptyScaffold() ){
		// add to parent removed list
		if( newContig->IsRepeat() )
			newScaffold->GetParent()->AddRemovedUnassignedNodes( newContig->GetOriginalRepeat(), newContigOri );
		else
			newScaffold->GetParent()->AddRemovedUnassignedNodes( newContig, newContigOri );
	}

	if( !newScaffold->GetParent()->IfEmptyScaffold() ){
		// check if the other orientation exists, remove and add to removed list
		int otherOri = (newContigOri+1) % 2;

		/*if( m_isSpecialGraph ){
			cerr<<newContig->GetName()<<"\t";
			if( newContigOri == PLUS )
				cerr<<"+\n";
			else
				cerr<<"-\n";
			
			if( otherOri == PLUS )
				cerr<<"other ori: +\n";
			else
				cerr<<"other ori: -\n";

			cerr<<"All unassigned nodes: \n";
			for( list< pair<Contig*, int> >::iterator iter = m_unassignedNodes->begin(); iter != m_unassignedNodes->end(); iter++ ){
				cerr<<(*iter).first->GetName()<<"\t";
				if( (*iter).second == PLUS )
					cerr<<"+\n";
				else 
					cerr<<"-\n";
			}
			cerr<<"\n";
		}
		*/

		/*cout<<newContig->GetName()<<endl;
		if( newContig->IsRepeat() )
			cout<<"is repeat\n";
		else
			cout<<"is not repeat\n";
		*/
		if( newContig->IsUnassignedNode( otherOri ) ){
			// add to current remove list
			//newScaffold->AddRemovedUnassignedNodes( newContig, otherOri );

			// remove from unassigned node list
			m_unassignedNodes->erase( newContig->GetUnassignedNodeIter( otherOri ) );
			newContig->ClearUnassignedNodeIter( otherOri );

			if( newContig->IsRepeat() ){
				newContig->GetOriginalRepeat()->ClearUnassignedNodeIter( otherOri );
			}
				
		}
	}

	//cerr<<"check connectivity\n";
	// check connectivity and update dangling edges and active region
	if( !CheckScaffold( newScaffold ) ){
		//fprintf( logFile, "Add fail\n" );
		//cerr<<"check scaffold fail\n";
		return false;
	}
	//cerr<<"finish checking\n";

	if( newScaffold->GetNumOfUnhappyEdges() > maxUnhappyEdge ){
		//fprintf( logFile, "Add fail\n" );
		//cerr<<"too many discordant edges: "<<newScaffold->GetNumOfUnhappyEdges()<<endl;
		return false;
	}

	// generate active region string and unhappy dangling edge string
	GenUDEString( newScaffold );
	GenARString( newScaffold );

	if( !m_activeRegion->empty() ){
		// calculate unassigned nodes
		//cerr<<"calcualte the unassigned nodes\n";
		//GenUnassigndNodes( newScaffold, false, false );

		/*cerr<<"Previous unassigned nodes:\n";
		for( list< pair<Contig*, int> >::iterator iter = m_unassignedNodes->begin(); iter != m_unassignedNodes->end(); iter++ ){
			cerr<<(*iter).first->GetName()<<"\t"<<(*iter).second<<endl;
		}
		cerr<<endl;
		*/
		
		
		// if need to recalculate
		// save the previous unassigned nodes
		//newScaffold->SavePreviousUnassignedNodesOrder( m_unassignedNodes );
		
		// recalculate the unassigned nodes
		//cerr<<"recalculate unassigned nodes.\n";
		GenUnassigndNodes( newScaffold, false, true );

		/*cerr<<"New unassigned nodes:\n";
		for( list< pair<Contig*, int> >::iterator iter = m_unassignedNodes->begin(); iter != m_unassignedNodes->end(); iter++ ){
			cerr<<(*iter).first->GetName()<<"\t"<<(*iter).second<<endl;
		}
		cerr<<endl;
		*/
		
	}
	else{
		// clear unassigned nodes
		//cerr<<"clear unassigned nodes\n";
		for( list< pair<Contig*, int> >::iterator iter = m_unassignedNodes->begin(); iter != m_unassignedNodes->end(); iter++ ){
			newScaffold->AddRemovedUnassignedNodes( (*iter).first, (*iter).second );
			(*iter).first->ClearUnassignedNodeIter( (*iter).second );
		}
		m_unassignedNodes->clear();
	}

	//cerr<<"end of adding contig"<<endl;

	return true;
}


// check edges of newly added contig
// ori is the orientation of contig c
void opera::CheckEdgesOfNewContig( PartialScaffold *p, Contig *c, list<PET*> *edges, int ori ){
	//cerr<<"check edges of new contig: "<<c->GetName()<<endl;
	for( list<PET*>::iterator iter = edges->begin(); iter != edges->end(); iter++ ){
		PET *edge = *iter;
		//cerr<<edge->GetOriString()<<endl;
		if( !(*iter)->isInSubgraph() )
			continue;

		if( (*iter)->GetUniqueEdge() != NULL || (*iter)->GetRepetitiveEdge() != NULL ){
			// do not consider the used edge
			//cerr<<"The edge is used: "<<(*iter)->GetOriString()<<endl;
			continue;
		}

		if( !c->IsRepeat() ){
			// c is unique contig
			if( edge->IsDE() ){
				// current edge is dangling edge
				if( !edge->IfUnhappy() ){
					// current edge was happy, remove from the list
					p->AddRemovedHappyDE( edge );
					m_happyDanglingEdges->erase( edge->GetHappyDEIter() );
					
					//if it is not happy any more, add the unhappy edge number
					if( !CheckDanglingEdgeHappiness( edge, c, ori ) ){
						p->AddUnhappyEdgeNumber( 1 );
						p->AddAddedUnhappyEdge( edge );
						//cerr<<"A: found an unhappy edge: "<<edge->GetOriString()<<endl;
					}
				}
				else{
					// current edge was not happy before, 
					// remove from unhappy dangling edges
					p->AddRmovedUnhappyDE( edge );
					m_unhappyDanglingEdges->erase( edge->GetUnhappyDEIter() );
				}
			}
			else{
				Contig *otherContig = edge->GetOtherContig( c );
				// it was not a dangling edge, check if orientation is correct
				if( !IsLeftEdge( edge, c, ori ) ){
					if( edge->GetDis() + edge->GetStd() * Configure::STD_TIMES >= -Configure::KMER ){
						// The edge extends to the right, it is happy, set as happy dangling edge
						p->AddAddedHappyDE( edge );
						m_happyDanglingEdges->push_front( edge );
						edge->SetHappyDEIter( m_happyDanglingEdges->begin() );

						
					}
				}
				else{
					// The edge extends to the left
					if( otherContig->IsRepeat() ){
						// FIXME: if c is connected to a repeat, check if it can be satisfied or not
						//cout<<edge->GetStartContig()->GetName()<<"\t"<<edge->GetEndContig()->GetName()<<endl;
						//cout<<"repeat is "<<otherContig->GetName()<<endl;
						int otherOri = edge->GetOrientationOfContig( otherContig );
						if( edge->GetPositionOfContig( otherContig ) == END )
							otherOri = GetOppositeOri( otherOri );
						
						// check and create repeat edge
						if( !CheckRepeatEdgeOfNewContig( otherContig, otherOri, c, edge ) ){
							// this edge is not happy 
							p->AddUnhappyEdgeNumber( 1 );
							p->AddAddedUnhappyEdge( edge );
							//cerr<<"B: found an unhappy edge: "<<edge->GetOriString()<<endl;
						}
						else{
							// FIXME: create repetitive edge
							//AddRepetitiveEdges( edge, otherContig );
						}
					}
					else{
						// if c is connected to a unique node, the edge is unhappy, set as unhappy dangling edge
						//cout<<"wrong orientation:"<<edge->GetOriString()<<endl;
						p->AddAddedUnhappyDE( edge, edge->GetOtherContig( c ) );
						m_unhappyDanglingEdges->push_front( edge );
						edge->SetUnhappyDEIter( m_unhappyDanglingEdges->begin() );
						//cerr<<"C: found an unhappy edge: "<<edge->GetOriString()<<endl;
					}
				}
			}
		}
		else{
			//cerr<<c->GetName()<<" is repeat\n";
			//cerr<<"edge is: "<<edge->GetOriString()<<endl;
			// c is repeat
			if( edge->IsDE() ){
				//cerr<<" it is dangling edge\n";
				// current edge is dangling edge
				if( !edge->IfUnhappy() ){
					//cerr<<"it was happy\n";
					// current edge was happy
					// FIXME: check if this edge can be satisfied by the newly added contig
					// fixed
					if( IsLeftEdge( edge, c->GetOriginalRepeat(), ori ) ){
						//cerr<<"it is satisfied\n";
						// if c has correct orientation, then the edge is satisfied
						p->AddRemovedHappyDE( edge );
						m_happyDanglingEdges->erase( edge->GetHappyDEIter() );

						// FIXME: create repeat edge
						// fixed
						//AddRepetitiveEdges( edge, c );
					}
					else{
						// if c does not have correct orientation, skip
						//cerr<<"it is not satisfied\n";
					}
				}
				else{
					// current edge was not happy before, 
					// remove from unhappy dangling edges
					p->AddRmovedUnhappyDE( edge );
					m_unhappyDanglingEdges->erase( edge->GetUnhappyDEIter() );
				}
			}
		}
	}
}

// check the satisfaction of repeated related edge of a new contig
bool opera::CheckRepeatEdgeOfNewContig( Contig *repeatContig, int repeatOri, Contig *addedContig, PET *edge ){
	int upperbound = edge->GetDis() + Configure::STD_TIMES * edge->GetStd();
	
	int dis = 0;
	int disWithGap = 0;
	SUBCONTIG *previousHead = NULL;
	for( list<Contig*>::reverse_iterator iter = m_activeRegion->rbegin(); iter != m_activeRegion->rend(); iter++ ){
		if( !(*iter)->IsRepeat() ){
			list<SUBCONTIG*> *subcontigs = (*iter)->GetSubcontigs();
			double length = 0;
			for( list<SUBCONTIG*>::reverse_iterator subIter = subcontigs->rbegin(); subIter != subcontigs->rend(); subIter++ ){
				if( (*subIter)->m_contigName == repeatContig->GetName() && (*subIter)->m_ori == repeatOri ){
					// this edge is satisfied
					// FIXME: create the unique edge for current edge
					// fixed
					//cerr<<"create unique edge in checkRepeatEdgeOfNewContig\n";
					CreateUniqueEdge( (*iter), (*iter)->GetOri(), addedContig, addedContig->GetOri(), length, disWithGap, edge );
					return true;
				}

				if( previousHead != NULL && subIter == subcontigs->rbegin() 
				    && (*subIter)->m_contigName == previousHead->m_contigName && (*subIter)->m_ori == previousHead->m_ori ){
					// ignore the tail of this contig if it is the same as the head of previous contig
					continue;
				}
				
				length += (*subIter)->m_length;
				
				dis += (*subIter)->m_length;
				disWithGap += (*subIter)->m_length;
				disWithGap += (*subIter)->m_gapSize;
				if( dis > upperbound + Configure::KMER )
					return false;
			}

			previousHead = *(subcontigs->begin());
		}
		else{
			// current contig is repeat, directly check
			if( (*iter)->GetName() == repeatContig->GetName() && (*iter)->GetOri() == repeatOri ){
				// FIXME: create repeat edge
				return true;
			}
			
			dis += (*iter)->GetLength();
			if( dis > upperbound + Configure::KMER )
				return false;
		}
	}

	return false;
}

// create the unique edge of an edge connecting with a repeat
void opera::CreateUniqueEdge( Contig *startContig, int startOri, Contig *endContig, int endOri, double disDelta, double disDeltaWithGap, PET *oriEdge ){
	//cerr<<"Create unique edge\n";
	PET *newEdge = new PET( startContig, startOri, endContig, endOri, oriEdge->GetDis() - disDelta, oriEdge->GetStd(), oriEdge->GetSize() );
	//cerr<<"set unique edge\n";
	oriEdge->SetUniqueEdge( newEdge );
	//cerr<<"set original edge\n";
	newEdge->SetOriginalEdge( oriEdge );
	//cerr<<"set repetitive edge as NULL\n";
	oriEdge->SetRepetitiveEdge( NULL );
	newEdge->SetDisWithGap( oriEdge->GetDisWithGap() - disDeltaWithGap );
	//cerr<<"successful\n";
}

// check dangling edge happiness
// newContig is the contig which is not in active region
// ori is the orientation of newContig
bool opera::CheckDanglingEdgeHappiness( PET *edge, Contig *newContig, int ori ){
	// check if edge is left edge( orientation is happy )
	if( !IsLeftEdge( edge, newContig, ori ) )
		return false;

	// check if distance is happy
		// first, find the node in active region
	Contig *activeNode = edge->GetOtherContig( newContig );
	//Contig *endContig = *m_activeRegion->rbegin();		// the last contig in active region

	list<Contig*>::reverse_iterator countIter = m_activeRegion->rbegin();
	int count = 0;
	while( *countIter != activeNode ){
		count++;
		countIter++;
	}

	int dis = (int) (newContig->GetStartPosition() - activeNode->GetStartPosition() - activeNode->GetLength() - count*Configure::KMER);
	if( dis - Configure::KMER > edge->GetDis() + Configure::STD_TIMES * edge->GetStd() ){
		return false;
	}
	return true;
}

// check if it is a left edge of current contig with certain orientation
bool opera::IsLeftEdge( PET *edge, Contig *newContig, int ori ){
	int petOri = edge->GetOrientationOfContig( newContig );
	int petPos = edge->GetPositionOfContig( newContig );
	if( ( petOri == ori && petPos == START ) || ( petOri != ori && petPos == END ) )
		return false;
	else
		return true;
}

// check connectivity and update scaffold active region & dangling edges
// return false if it is not connected
bool opera::CheckScaffold( PartialScaffold *p ){
	Contig *endContig = *m_activeRegion->rbegin();		// the last contig in active region

	//cerr<<"Check dangling edge\n";
	m_repeatTraceBack = NULL;
	
	// if the added contig is not repeat, check each dangling edge with repeat to see if it can be satisfied by the newly added node
	// FIXME: if an edge with repeat is satisfied, create a unique edge
	// fixed
	if( !endContig->IsRepeat() ){
		p->GetAddedContigInAR()->SplitSuperContig( p->GetParent()->GetAddedContigInAR() );
		int removedDis = p->GetAddedContigInAR()->GetRemovedDis();
	
		// remove one occurrence if a repeat occured twice consecutively
		endContig->SetStartPosition( endContig->GetStartPosition() - removedDis, 0 );

		// FIXME: if two occurrence of the same repeat are adjacent, remove one and reduce the distance
		// fixed
		list<PET*>::iterator edgeIter = m_happyDanglingEdges->begin();
		while( edgeIter != m_happyDanglingEdges->end() ){
			PET *edge = *edgeIter;
			// determine if this edge connects with repeat
			// get the contig in active region
			Contig *ARContig = edge->FindARContig();
			if( ARContig->IsRepeat() ){
				continue;
			}

			Contig *otherContig = edge->GetOtherContig( ARContig );
			if( !otherContig->IsRepeat() ){
				edgeIter++;
				continue;
			}
			
			// determine if the edge can be satisfied or not
			//cerr<<p->GetAddedContigInAR()->GetName()<<endl;
			//cerr<<edge->GetOriString()<<endl;
			Contig *uniqueARContig = NULL;
			int uniqueDisDelta;
			int uniqueDisDeltaWithGap;
			if( !p->GetAddedContigInAR()->IsDanglingEdge( uniqueARContig, uniqueDisDelta, uniqueDisDeltaWithGap, edge ) ){
				//
				// FIXME: Create a unique edge for current edge
				// fixed: 
				//cerr<<"create a unqiue edge for current edge\n";
				//if( uniqueARContig == NULL )
				//	cerr<<"ar contig is null\n";
				//cerr<<uniqueARContig->GetName()<<endl;
				CreateUniqueEdge( uniqueARContig, uniqueARContig->GetOri(), 
						  p->GetAddedContigInAR(), p->GetAddedContigInAR()->GetOri(), uniqueDisDelta, uniqueDisDeltaWithGap, edge );
				edgeIter = m_happyDanglingEdges->erase( edgeIter );
				edge->SetDE( false );
				p->AddRemovedHappyDE( edge );
				//cerr<<"it is not a dangling edge any more\n";
			}
			else{
				edgeIter++;
				//cerr<<"It is still a dangling edge\n";
			}
		}
	}
	//////////////////
	//cerr<<p->GetNumOfUnhappyEdges()<<endl;

	//cerr<<"check distance\n";
	// check each dangling edge if the distance could be satisfied
	//////FIXME: if two occurrence of the same repeat are adjacent, remove one and reduce the distance
	// FIXED: remove the corresponding distance
	list<Contig*>::iterator arIter = m_activeRegion->begin();
	//cerr<<"Active region:\n";
	while( arIter != m_activeRegion->end() ){
		//cerr<<(*arIter)->GetName()<<endl;
		arIter++;
	}
	//cerr<<"#######"<<endl;

	//cerr<<"Check dangling edge\n";
	list<PET*>::iterator iter = m_happyDanglingEdges->begin();
	while( iter != m_happyDanglingEdges->end() ){
		PET *edge = *iter;
		//cerr<<edge->GetOriString()<<endl;
		//if( (*iter)->IsDE() && !(*iter)->IfUnhappy() && (*iter)->isInSubgraph() )
		//	{}
		//else
		//	cerr<<" This edge should be happy!\n";
		/*if( (*iter)->IsDE() )
			cerr<<"it is de\n";
		if( !(*iter)->IfUnhappy() )
			cerr<<"it is happy\n";
		if( (*iter)->isInSubgraph() ) 
		cerr<<"it is in subgraph\n";*/
		

		// get the contig in active region
		Contig *ARContig = edge->FindARContig();
		//cerr<<"AR contig is "<<ARContig->GetName()<<endl;
		list<Contig*>::reverse_iterator countIter = m_activeRegion->rbegin();
		int num = 0;	// record the number of contigs lying in this edge
		while( *countIter != ARContig ){
			num++;
			countIter++;
		}
		//cerr<<"found ar contig\n";

		//int dis = (int) (endContig->GetStartPosition() + endContig->GetLength() - endContig->GetRemovedDis()
		//		 - ARContig->GetStartPosition() - ARContig->GetLength() - ARContig->GetRemovedDis()- 2) - num * Configure::KMER;
		int dis = (int) (endContig->GetStartPosition() + endContig->GetLength() - endContig->GetRemovedDis()
				 - ARContig->GetStartPosition() - ARContig->GetLength() - ARContig->GetRemovedDis()- 2) + Configure::KMER;
		/*if( endContig->GetName() == "462051" && ARContig->GetName() == "499283" ){
			cerr<<"find Edge: \n\tdis: "<<dis<<endl;
			cerr<<"\t"<<endContig->GetStartPosition()<<"\t"<<endContig->GetLength()<<"\t"<<endContig->GetRemovedDis()<<endl;
			cerr<<"\t"<<ARContig->GetStartPosition()<<"\t"<<ARContig->GetLength()<<"\t"<<ARContig->GetRemovedDis()<<endl;
			cerr<<"******************************************\n";
		}
		*/

		//cerr<<"left: "<<dis - 2*Configure::KMER<<endl;
		//cerr<<"right: "<<edge->GetDis() + Configure::STD_TIMES * edge->GetStd()<<endl;
		if( dis - 2*Configure::KMER >= edge->GetDis() + Configure::STD_TIMES * edge->GetStd() ){
			// change to unhappy dangling edge
			//cerr<<"set to unhappy dangling edge\n";
			iter = m_happyDanglingEdges->erase( edge->GetHappyDEIter() );
			p->AddRemovedHappyDE( edge );
			//cout<<"distance not happy:"<<edge->GetOriString()<<endl;
			// only set the edges with unique nodes as unhappy dangling edges
			if( !edge->GetOtherContig( ARContig )->IsRepeat()){
				m_unhappyDanglingEdges->push_front( edge );
				edge->SetUnhappyDEIter( m_unhappyDanglingEdges->begin() );
				p->AddAddedUnhappyDE( edge, edge->GetOtherContig( ARContig ) );
			}
			else{
				// for repeat, only set this edge as unhappy edge
				p->AddUnhappyEdgeNumber( 1 );
				p->AddAddedUnhappyEdge( edge );
			}
			//cerr<<edge->GetOriString()<<endl;
			//cerr<<"dis: "<<dis<<"\n";
		}
		else
			iter++;
	}
	//cerr<<p->GetNumOfUnhappyEdges()<<endl;

	// update the active region, remove the head nodes which do not have happy dangling edges
	//cerr<<"Update active region\n";
	// without repeat
	/*list<Contig*>::iterator contigIter = m_activeRegion->begin();
	while( contigIter != m_activeRegion->end() ){
		Contig *currentContig = *contigIter;
		if( !currentContig->IfHasHappyDE( currentContig->GetOri() ) ){
			// remove from active region
			//cerr<<"Remove nodes from active region: "<<currentContig->GetName()<<endl;
			contigIter = m_activeRegion->erase( contigIter );
			currentContig->SetIfInAR( false );
			if( currentContig != p->GetAddedContigInAR() )
				p->AddRemovedContigInAR( currentContig );
			else
				p->RemoveOneFromActiveRegion();		// it is the newly added contig, 
													// only minus 1 from number of active region

			// check if currentContig connects to active region
			//if( !CheckActiveRegionConnectivity( currentContig ) )
			//	return false;
			
			//cout<<"remove node from ar successfully: "<<currentContig->GetName()<<endl;
		}
		else
			break;
	}
	*/

	//cerr<<"remove contigs without edges\n";
	// with repeat
	list<Contig*>::iterator contigIter = m_activeRegion->begin();
	while( contigIter != m_activeRegion->end() ){
		Contig *currentContig = *contigIter;
		if( !currentContig->IsRepeat() && !currentContig->IfHasHappyDE( currentContig->GetOri() ) ){
			// if it is a unique contig without any happy dangling edges
			// remove from active region
			//cerr<<"Remove nodes from active region: "<<currentContig->GetName()<<endl;
			contigIter = m_activeRegion->erase( contigIter );
			currentContig->SetIfInAR( false );
			if( currentContig != p->GetAddedContigInAR() )
				p->AddRemovedContigInAR( currentContig );
			else
				p->RemoveOneFromActiveRegion();		// it is the newly added contig, 
													// only minus 1 from number of active region

			// check if currentContig connects to active region
			//if( !CheckActiveRegionConnectivity( currentContig ) )
			//	return false;
			
			//cout<<"remove node from ar successfully: "<<currentContig->GetName()<<endl;
		}
		else
			break;
	}
	//cerr<<p->GetNumOfUnhappyEdges()<<endl;

	// shrink the active region if it is longer than library threshold
	if( m_activeRegion->empty() ){
		return true;
	}

        //cerr<<"shrink active region\n";
	list<Contig*>::reverse_iterator rIter = m_activeRegion->rbegin();
	double length = 0;
	while( rIter != m_activeRegion->rend() ){
		Contig *currentContig = *rIter;
		length += currentContig->GetLength() - Configure::KMER - currentContig->GetRemovedDis();
		if( length - Configure::KMER > Configure::UPPERBOUND ){ //Configure::LIB_MEAN + Configure::STD_TIMES * Configure::LIB_STD ){
			// have to remove all contigs before current contig
			rIter++;
			break;
		}
		rIter++;
	}

	// FIXME: if all the contigs in the active region do not have dangling edges, remove everything
	list<Contig*>::reverse_iterator rIterOther = m_activeRegion->rbegin();
	while( rIterOther != rIter ){
		Contig *currentContig = *rIterOther;
		if( !currentContig->IsRepeat() && currentContig->IfHasHappyDE( currentContig->GetOri() )){
			// if current contig is unique and has happy dangling edges, don't need to empty the active region
			break;
		}
		rIterOther++;
	}

	if( rIterOther == rIter ){
		// need to empty the active region
		rIter = m_activeRegion->rbegin();
	}

	if( rIter == m_activeRegion->rend() )
		return true;

	Contig *lastRemovedContig = *rIter;

	

	// return to the starting contig
	//rIter--;

	//cerr<<"modify active region edge\n";
	contigIter = m_activeRegion->begin();
	int posInAR = m_activeRegion->size() + 1;
	list<Contig*> erasedARContigs;   // record all the contigs erased from the active region
	int tempNum = 0;
	int ARSize = m_activeRegion->size();
	while( contigIter != m_activeRegion->end() ){ //&& contigIter != rIter ){
		tempNum++;
		//cerr<<"deleting "<<tempNum<<" contig\n";
 		Contig *currentContig = *contigIter;
		posInAR--;

		//cerr<<"deleting "<<currentContig->GetName()<<endl;
		
		// remove from active region
		contigIter = m_activeRegion->erase( contigIter );
		currentContig->SetIfInAR( false );           ///////////FIXME: repeat has multiple occurrences
		                                             ///////////FIXED: by traversing active region later and labeling all contigs
		erasedARContigs.push_back( currentContig );

		if( currentContig->IsRepeat() ){
			// FIXME: check if this repeat can be removed without introducing any discordant edges or not
			// FIXED: added the code
			//cerr<<"repeat\n";
			bool needToRemove = BringDiscordantEdgeByRemoving( currentContig, currentContig->GetOri(), erasedARContigs, p, ARSize);
			if( needToRemove ){
				m_repeatTraceBack = currentContig;
				return false;
			}
			//cerr<<"finish repeat\n";
			//////////////////////
		}
		else{
			//cerr<<"unique contig\n";
			// remove from active region
			//contigIter = m_activeRegion->erase( contigIter );
			//currentContig->SetIfInAR( false );           ///////////FIXME: repeat has multiple occurrences
			                                             ///////////FIXED: by traversing active region later and labeling all contigs
			
			// set all happy dangling edges to unhappy dangling edges
			list<PET*>::iterator edgeIter = currentContig->GetLeftEdges()->begin();
			while( edgeIter != currentContig->GetLeftEdges()->end() ){
				if( (*edgeIter)->IsDE() && !(*edgeIter)->IfUnhappy() ){
					if( !(*edgeIter)->GetOtherContig( currentContig )->IsRepeat() ){
						// change to unhappy dangling edge
						m_happyDanglingEdges->erase( (*edgeIter)->GetHappyDEIter() );
						m_unhappyDanglingEdges->push_front( *edgeIter );
						(*edgeIter)->SetUnhappyDEIter( m_unhappyDanglingEdges->begin() );
						p->AddRemovedHappyDE( *edgeIter );
						p->AddAddedUnhappyDE( *edgeIter, (*edgeIter)->GetOtherContig( currentContig ) );
					}
					else{
						// change to unhappy edge
						m_happyDanglingEdges->erase( (*edgeIter)->GetHappyDEIter() );
						p->AddRemovedHappyDE( *edgeIter );
						p->AddUnhappyEdgeNumber( 1 );
						p->AddAddedUnhappyEdge( *edgeIter );
					}
				}
				edgeIter++;
			}
			edgeIter = currentContig->GetRightEdges()->begin();
			while( edgeIter != currentContig->GetRightEdges()->end() ){
				if( (*edgeIter)->IsDE() && !(*edgeIter)->IfUnhappy() ){
					if( !(*edgeIter)->GetOtherContig( currentContig )->IsRepeat() ){
						// change to unhappy dangling edge
						m_happyDanglingEdges->erase( (*edgeIter)->GetHappyDEIter() );
						m_unhappyDanglingEdges->push_front( *edgeIter );
						(*edgeIter)->SetUnhappyDEIter( m_unhappyDanglingEdges->begin() );
						p->AddRemovedHappyDE( *edgeIter );
						p->AddAddedUnhappyDE( *edgeIter, (*edgeIter)->GetOtherContig( currentContig ) );
					}
					else{
						// change to unhappy edge
						m_happyDanglingEdges->erase( (*edgeIter)->GetHappyDEIter() );
						p->AddRemovedHappyDE( *edgeIter );
						p->AddUnhappyEdgeNumber( 1 );
						p->AddAddedUnhappyEdge( *edgeIter );
					}
				}
				edgeIter++;
			}
			//cerr<<"finish unique\n";
		}
			
		if( currentContig != p->GetAddedContigInAR() )
			p->AddRemovedContigInAR( currentContig );
		else
			p->RemoveOneFromActiveRegion();		// it is the newly added contig, 
		// only minus 1 from number of active region

		// check if currentContig connects to active region
		//if( !CheckActiveRegionConnectivity( currentContig ) )
		//	return false;

		/*if( currentContig->GetName() == "1455" ){
			// print the edge
			for( list<PET*>::iterator tempEdgeIter = currentContig->GetRightEdges()->begin(); 
			     tempEdgeIter != currentContig->GetRightEdges()->end(); tempEdgeIter++ ){
				cerr<<(*tempEdgeIter)->GetOriString()<<endl;
				cerr<<(*tempEdgeIter)<<endl;
			}
			}*/

		if( currentContig == lastRemovedContig ){
			// this is the last removed contig, do not continue to traverse
			break;
		}
	}
	//cerr<<p->GetNumOfUnhappyEdges()<<endl;
	
	// set all contigs in active region as being in AR
	/*contigIter = m_activeRegion->begin();
	while( contigIter != m_activeRegion->end() ){
		Contig *currentContig = *contigIter;
		currentContig->SetIfInAR( true );
	}
	*/
	
	//cerr<<"finish\n";
	return true;
}



// check if contig c could connect to active region
// return false if not connected
bool opera::CheckActiveRegionConnectivity( Contig *c ){
	if( m_activeRegion->empty() ){
		//cout<<"return from empty"<<endl;
		return true;
	}

	set<Contig*> *visitedContigs = new set<Contig*>;
	list<Contig*> *possibleContigs = new list<Contig*>;

	possibleContigs->push_back( c );
	visitedContigs->insert( c );

	while( !possibleContigs->empty() ){
		Contig *currentContig = *possibleContigs->begin();
		possibleContigs->pop_front();

		// check left edges
		list<PET*> *edges = currentContig->GetLeftEdges();
		for( list<PET*>::iterator iter = edges->begin(); iter != edges->end(); iter++ ){
			if( !(*iter)->isInSubgraph() || (*iter)->IfUnhappy() )
				continue;

			Contig *otherContig = (*iter)->GetOtherContig( currentContig );
			if( otherContig->IfInAR() ){
				//cout<<"The contig "<<c->GetName()<<" which is connected with in active region: "<<otherContig->GetName()<<endl;
				delete visitedContigs;
				delete possibleContigs;
				return true;
			}

			if( visitedContigs->find( otherContig ) == visitedContigs->end() ){
				visitedContigs->insert( otherContig );
				possibleContigs->push_back( otherContig );
			}
		}

		// check right edges
		edges = currentContig->GetRightEdges();
		for( list<PET*>::iterator iter = edges->begin(); iter != edges->end(); iter++ ){
			if( !(*iter)->isInSubgraph() || (*iter)->IfUnhappy() )
				continue;

			Contig *otherContig = (*iter)->GetOtherContig( currentContig );
			if( otherContig->IfInAR() ){
				//cout<<"The contig "<<c->GetName()<<" which is connected with in active region: "<<otherContig->GetName()<<endl;
				delete visitedContigs;
				delete possibleContigs;
				return true;
			}

			if( visitedContigs->find( otherContig ) == visitedContigs->end() ){
				visitedContigs->insert( otherContig );
				possibleContigs->push_back( otherContig );
			}
		}
	}

	delete visitedContigs;
	delete possibleContigs;
	return false;
}

// recover unassigned nodes
void opera::RecoverUnassignedNodes( list< pair<Contig*, int> > *nodes ){
	for( list< pair<Contig*, int> >::iterator iter = m_unassignedNodes->begin(); iter != m_unassignedNodes->end(); iter++ ){
		(*iter).first->ClearUnassignedNodeIter( (*iter).second );
	}

	m_unassignedNodes->clear();
	
	for( list< pair<Contig*, int> >::reverse_iterator iter = nodes->rbegin(); iter != nodes->rend(); iter++ ){
		m_unassignedNodes->push_front( *iter );
		(*iter).first->SetUnassignedNodeIter( (*iter).second, m_unassignedNodes->begin() );
	}
}


// trace back to parent partial scaffold
// change p points to its parent
void opera::TraceBack( PartialScaffold *&p ){
	//cerr<<"trace back\n";
	bool continueToTraceBack = false;

	// if the unassigned nodes have been recalculated, recover the unassigned nodes
	//if( p->IfRecalculatdUnassignedNodes() ){
		RecoverUnassignedNodes( p->GetPreviousUnassignedNodesOrder() );
		//}
	
	// recover active region
	Contig *addedContig = p->GetAddedContigInAR();
	if( m_repeatTraceBack != NULL ){
		if( addedContig != m_repeatTraceBack )
			continueToTraceBack = true;
		else
			m_repeatTraceBack = NULL;   // reset the repeat tracing back to
	}

	/*if( addedContig != NULL && addedContig->IfInAR() ){
		addedContig->SetIfInAR( false );
		m_activeRegion->pop_back();
		
		// if it is a repeat, delete this copy
		if( addedContig->IsRepeat() )
			delete addedContig;
		else{
			// if it is a unique node, clear the status of all edges
			addedContig->ClearStatusOfEdges( addedContig->GetLeftEdges() );
			addedContig->ClearStatusOfEdges( addedContig->GetRightEdges() );
		}
	}
	*/

	//*************modified
	if( addedContig != NULL ){
		if( addedContig->IfInAR() ){
			addedContig->SetIfInAR( false );
			m_activeRegion->pop_back();
		}
		
		// if it is a repeat, delete this copy
		if( addedContig->IsRepeat() ){
			//addedContig->DeleteEdgesForRepeat( addedContig->GetLeftEdges() );
			//addedContig->DeleteEdgesForRepeat( addedContig->GetRightEdges() );
			delete addedContig;
		}
		else{
			// if it is a unique node, clear the status of all edges
			addedContig->ClearStatusOfEdges( addedContig->GetLeftEdges() );
			addedContig->ClearStatusOfEdges( addedContig->GetRightEdges() );
		}
	}
	//****************************

	list<Contig*> *removedContigs = p->GetRemovedContigsInAR();
	for( list<Contig*>::reverse_iterator iter = removedContigs->rbegin(); 
		iter != removedContigs->rend(); iter++ ){
			m_activeRegion->push_front( *iter );
			(*iter)->SetIfInAR( true );
			
			if( (*iter)->IsRepeat() ){
				// this contig is repeat, remove all its edges
				(*iter)->DeleteEdgesForRepeat( (*iter)->GetLeftEdges() );
				(*iter)->DeleteEdgesForRepeat( (*iter)->GetRightEdges() );
			}
	}

	

	// recover added happy and unhappy dangling edges
	list<PET*> *addedPets = p->GetAddedHappyDE();
	list<PET*>::iterator iter;
	for( iter = addedPets->begin(); iter != addedPets->end(); iter++ ){
		PET *edge = *iter;
		m_happyDanglingEdges->erase( edge->GetHappyDEIter() );
		edge->SetDE( false );
	}

	addedPets = p->GetAddedUnhappyDE();
	for( iter = addedPets->begin(); iter != addedPets->end(); iter++ ){
		PET *edge = *iter;
		m_unhappyDanglingEdges->erase( edge->GetUnhappyDEIter() );
		edge->SetDE( false );
		edge->SetHappy();
	}

	// recover removed happy and unhappy dangling edges
	list<PET*> *removedPets = p->GetRemovedHappyDE();
	for( iter = removedPets->begin(); iter != removedPets->end(); iter++ ){
		PET *edge = *iter;
		//cerr<<"recover removed happy dangling edges: "<<edge->GetOriString()<<endl;
		m_happyDanglingEdges->push_front( edge );
		edge->SetHappyDEIter( m_happyDanglingEdges->begin() );
		edge->SetDE( true );
	}

	removedPets = p->GetRemovedUnhappyDE();
	for( iter = removedPets->begin(); iter != removedPets->end(); iter++ ){
		PET *edge = *iter;
		m_unhappyDanglingEdges->push_front( edge );
		edge->SetUnhappyDEIter( m_unhappyDanglingEdges->begin() );
		edge->SetDE( true );
		edge->SetUnhappy();
	}
	

	// recover unhappy edges
	addedPets = p->GetAddedUnhappyEdges();
	for( iter = addedPets->begin(); iter != addedPets->end(); iter++ ){
		(*iter)->SetHappy();
	}

	// recover unassigned node sets
	/*list< pair<Contig*, int> >::iterator nodeIter;
	list< pair<Contig*, int> > *removedNodes = p->GetRemovedUnassignedNodes();
	for( nodeIter = removedNodes->begin(); nodeIter != removedNodes->end(); nodeIter++ ){
		pair<Contig*, int> newPair = *nodeIter;
		//m_unassignedNodes->push_front( newPair );
		//newPair.first->SetUnassignedNodeIter( newPair.second, m_unassignedNodes->begin() );
		InsertContigIntoUnassignedNodeSet( newPair );
	}

	list< pair<Contig*, int> > *addedNodes = p->GetAddedUnassignedNodes();
	for( nodeIter = addedNodes->begin(); nodeIter != addedNodes->end(); nodeIter++ ){
		pair<Contig*, int> newPair = *nodeIter;
		m_unassignedNodes->erase( newPair.first->GetUnassignedNodeIter( newPair.second ) );
		newPair.first->ClearUnassignedNodeIter( newPair.second );
	}
	*/

	PartialScaffold *tempSca = p;
	p = p->GetParent();
	delete tempSca;

	if( continueToTraceBack )
		TraceBack( p );
}

// generate results of scaffold
// save the scaffolds string
// and calculate the distance of each contigs in scaffolds
void opera::GenerateResults( PartialScaffold *p, ScaffoldResult **&results, int num ){
	//fprintf( logFile, "start generating results\n" );
	//fflush( logFile );

	// calculate gap
	PartialScaffold *tempPS = p;
	int scaffoldID = 0;
	int contigOrderID = 0;
	int numberOfGaps = 0;
	vector<Contig*> *contigInScaffold = new vector<Contig*>;		// save the contigs in one scaffold
	while( tempPS->GetParent() != NULL ){
		//if( tempPS->GetNumOfContigsInAR() == 0 && tempPS->GetNumOfHappyDE() == 0 ){
		if( tempPS->IfEmptyScaffold() ){
			// it is a new scaffold
			//if( tempPS->GetNumOfUnhappyDE() != 0 || tempPS->IfMannualEmptyScaffold() ){
				// calculate the gap sizes
				if( contigInScaffold->size() > 1 )
					CalculateGap( contigInScaffold, numberOfGaps - 1 );
				else{
#ifdef DEBUG
					cout<<"singleton scaffold"<<endl;
#endif
					contigInScaffold->at( 0 )->SetGap( 0 );
				}
				scaffoldID++;
			//}
			
			// start a new scaffold
			if( tempPS->IfMannualEmptyScaffold() )
				tempPS = tempPS->GetParent();

			contigOrderID = 0;
			numberOfGaps = 0;
			contigInScaffold->clear();
		}

		Contig *addedContig = tempPS->GetAddedContigInAR();
		addedContig->SetScaffoldPos( contigOrderID++ );
		addedContig->SetScaffoldID( scaffoldID );
		contigInScaffold->push_back( addedContig );
		numberOfGaps++;
		
		tempPS = tempPS->GetParent();
	}

	// calculate gap for last scaffold
	if( contigInScaffold->size() > 1 )
		CalculateGap( contigInScaffold, numberOfGaps - 1 );
	else{
#ifdef DEBUG
		cout<<"singleton scaffold"<<endl;
#endif
		contigInScaffold->at( 0 )->SetGap( 0 );
	}

	delete contigInScaffold;

	//fprintf( logFile, "finish calculating gaps\n" );
	//fflush( logFile );

	// generate scaffold
	//PartialScaffold *tail = p;
	string scaffoldResult = "";
	double scaffoldLength = 0;
	double scaffoldLengthWithGap = 0;
	double scaffoldCov = 0;		// the average coverage of scaffolds
	scaffoldID = 0;
	double endScaffoldPos = p->GetAddedContigInAR()->GetStartPosition()
					+ p->GetAddedContigInAR()->GetLength() - 1;
	double endScaffoldPosWithGap = p->GetAddedContigInAR()->GetStartPositionWithGap() 
		+ p->GetAddedContigInAR()->GetLengthWithGap() - 1;
	// cout<<"start to generate scaffold results"<<endl;
	bool firstContig = true;
	while( p->GetParent() != NULL ){
		// save unhappy edges information
		list<PET*> *unhappyEdges = p->GetAddedUnhappyEdges();
		unhappyEdges->insert( unhappyEdges->begin(), p->GetAddedUnhappyDE()->begin(),
			p->GetAddedUnhappyDE()->end() );
		for( list<PET*>::iterator iter = unhappyEdges->begin(); iter != unhappyEdges->end(); iter++ ){
			m_unhappyEdgeString->push_back( (*iter)->GetOriString() );
		}
		//if( p->GetNumOfContigsInAR() == 0 && p->GetNumOfHappyDE() == 0 ){
		if( p->IfEmptyScaffold() ){
			// it is a new scaffold
			//if( p->GetNumOfUnhappyDE() != 0  || p->IfMannualEmptyScaffold() ){
				// calculate the gap sizes
				// not the end of solution, save previous results
#ifdef SPLIT
				fprintf( logFile, "new scaffold result %d\n", scaffoldID );
				if( scaffoldID >= num )
					fprintf( logFile, "Error: scaffoldID %d > total number of scaffolds (%d)\n", scaffoldID, num );
#endif

				results[ scaffoldID ] = new ScaffoldResult();
				results[ scaffoldID ]->SetScaffoldString( scaffoldResult );
				results[ scaffoldID ]->SetLength( scaffoldLength );
				results[ scaffoldID ]->SetLengthWithGap( scaffoldLengthWithGap );
				scaffoldCov = scaffoldCov / scaffoldLength;
				results[ scaffoldID ]->SetCov( scaffoldCov );
				scaffoldID++;
			//}
			
			// start a new scaffold
			scaffoldResult = "";
			scaffoldLength = 0;
			scaffoldLengthWithGap = 0;
			scaffoldCov = 0;
			if( p->IfMannualEmptyScaffold() )
				p = p->GetParent();
			firstContig = true;

			endScaffoldPos = p->GetAddedContigInAR()->GetStartPosition()
					+ p->GetAddedContigInAR()->GetLength() - 1;
			endScaffoldPosWithGap = p->GetAddedContigInAR()->GetStartPositionWithGap() 
				+ p->GetAddedContigInAR()->GetLengthWithGap() - 1;
		}

		// generate the scaffold string
		Contig *addedContig = p->GetAddedContigInAR();
		if( !firstContig ){
			//scaffoldResult = "\n" + scaffoldResult;
		}
		else
			firstContig = false;

		// remove the adjacent same repeat
		//scaffoldResult = addedContig->GenScaffoldString( addedContig->GetOri() ) 
		//	+ itos( addedContig->GetGap() ) + scaffoldResult;
		cerr<<"add string to scaffold\n";
		scaffoldResult = AddScaffoldString( addedContig->GenScaffoldString( addedContig->GetOri() ) + itos( addedContig->GetGap() ), scaffoldResult );
		scaffoldLength += addedContig->GetLength();
		scaffoldLengthWithGap += addedContig->GetLengthWithGap() + addedContig->GetGap();
		scaffoldCov += addedContig->GetCov() * addedContig->GetLength();

		// calculate the left and right distance of contig
		addedContig->SetScaffoldID( scaffoldID );
		double leftDis = addedContig->GetStartPosition() - 1;
		double rightDis = endScaffoldPos - addedContig->GetStartPosition() - addedContig->GetLength() + 1;
		double leftDisWithGap = addedContig->GetStartPositionWithGap() - 1;
		double rightDisWithGap = endScaffoldPosWithGap - addedContig->GetStartPositionWithGap() - addedContig->GetLengthWithGap() + 1;
		//cout<<"Left: "<<leftDis<<"\t"<<leftDisWithGap<<endl;
		//cout<<"Right: "<<rightDis<<"\t"<<rightDisWithGap<<endl;
#ifdef CHECK
		if( leftDis < 0 ){
			cout<<addedContig->GetScaffoldString()<<endl;
			cout<<"leftDis < 0: "<<leftDis<<endl<<endl;
		}
		
		if( rightDis < 0 ){
			cout<<addedContig->GetScaffoldString()<<endl;
			cout<<"end scaffold pos: "<<endScaffoldPos<<endl;
			cout<<"rightDis < 0: "<<rightDis<<endl<<endl;
		}
#endif
		addedContig->SetLeftDis( leftDis );
		addedContig->SetRightDis( rightDis );
		addedContig->SetLeftDisWithGap( leftDisWithGap );
		addedContig->SetRightDisWithGap( rightDisWithGap );

		p = p->GetParent();
	}

	// save the last results
#ifdef SPLIT
	fprintf( logFile, "new scaffold result %d\n", scaffoldID );
#endif
	if( scaffoldID != num - 1 ){
		cout<<"ERROR: scaffoldID "<<scaffoldID<<" is not equal to total number of scaffolds "<<num<<endl;
		exit(0);
	}
	results[ scaffoldID ] = new ScaffoldResult();
	results[ scaffoldID ]->SetScaffoldString( scaffoldResult );
	results[ scaffoldID ]->SetLength( scaffoldLength );
	results[ scaffoldID ]->SetLengthWithGap( scaffoldLengthWithGap );
	scaffoldCov = scaffoldCov / scaffoldLength;
	results[ scaffoldID ]->SetCov( scaffoldCov );
}

// remove the labeled repeats
void opera::RemoveRepeats( vector<Contig*> *scafs ){
	// FIXME: after removal, modify the position in scaffold
	// scan the scaffolds
	//cerr<<"filter all the removed repeats\n";
	int index = 0;
	vector<Contig*>::iterator contigIter = scafs->begin();
	int pos = 0;
	while( contigIter != scafs->end() ){
		//cerr<<(*contigIter)->GetName()<<"\t"<<(*contigIter)->GetScaffoldPos()<<endl;
		//for( vector<Contig*>::iterator contigIter = scafs->begin(); contigIter != scafs->end(); contigIter++ ){
		// if current repeat need to be removed
		if( (*contigIter)->NeedRemove() ){
			//cerr<<"find one contig need to remove: "<<(*contigIter)->GetName()<<endl;
			Contig *repeat = *contigIter;
			int repeatOri = repeat->GetOri();

			// locate the previous and the next occurrences
			Contig *previousOccu = NULL;
			for( int i = index - 1; i >= 0; i-- ){
				if( scafs->at( i )->GetName() == repeat->GetName() && scafs->at( i )->GetOri() == repeatOri ){
					previousOccu = scafs->at( i );
					break;
				}
			}

			Contig *nextOccu = NULL;
			for( int i = index + 1; i < (int) scafs->size(); i++ ){
				if( scafs->at( i )->GetName() == repeat->GetName() && scafs->at( i )->GetOri() == repeatOri ){
					nextOccu = scafs->at( i );
					break;
				}
			}

			// assign edges to those two repeats
			if( previousOccu != NULL ){
				if( repeatOri == PLUS ){
					// assign the right edges
					AssignRepetitiveEdges( repeat->GetRightEdges(), previousOccu );
				}
				else{
					// assign the left edges
					AssignRepetitiveEdges( repeat->GetLeftEdges(), previousOccu );
				}
			}

			if( nextOccu != NULL ){
				if( repeatOri == PLUS ){
					// assign the left edges
					AssignRepetitiveEdges( repeat->GetLeftEdges(), nextOccu );
				}
				else{
					// assign the right edges
					AssignRepetitiveEdges( repeat->GetRightEdges(), nextOccu );
				}
			}

			// remove current repeat from the vector
			contigIter = scafs->erase( contigIter );
		}
		else{
			(*contigIter)->SetScaffoldPos( pos );
			pos++;
			contigIter++;
			index++;
		}
	}

	
}

// assign repetitive edges
void opera::AssignRepetitiveEdges( list<PET*> *edges, Contig *repeat ){
	//cerr<<"assign repetitive edges\n";
	for( list<PET*>::iterator iter = edges->begin(); iter != edges->end(); iter++ ){
		//int repeatPos;
		Contig *originalRepeat;
		if( (*iter)->GetStartContig()->GetName() == repeat->GetName() ){
			//repeatPos = PLUS;
			originalRepeat = (*iter)->GetStartContig();
		}
		else{
			//repeatPos = MINUS;
			originalRepeat = (*iter)->GetEndContig();
		}

		int repeatOri = (*iter)->GetOrientationOfContig( originalRepeat );
		Contig *otherContig = (*iter)->GetOtherContig( originalRepeat );
		int otherPos = (*iter)->GetPositionOfContig( otherContig );
		int otherOri = (*iter)->GetOrientationOfContig( otherContig );

		PET *newEdge;
		if( otherPos == START ){
			newEdge = new PET( otherContig, otherOri, repeat, repeatOri, (*iter)->GetDis(), (*iter)->GetStd(), (*iter)->GetSize() );
		}
		else{
			newEdge = new PET( repeat, repeatOri, otherContig, otherOri, (*iter)->GetDis(), (*iter)->GetStd(), (*iter)->GetSize() );
		}

		(*iter)->GetOriginalEdge()->SetRepetitiveEdge( newEdge );
		newEdge->SetOriginalEdge( (*iter)->GetOriginalEdge() );
		repeat->AddEdge( newEdge );
	}
}

void opera::GenerateResults( vector<vector<Contig*>*> *scafs, ScaffoldResult **&results ){
	// calculate gap
	for( int i = 0; i < (int) scafs->size(); i++ ){
		if( scafs->at( i )->size() == 1 ){
			//cout<<"singleton\n";
			scafs->at( i )->at( 0 )->SetGap( 0 );
		}
		else{
			/*cerr<<"non singleton\n";
			cerr<<"total contigs: "<<scafs->at( i )->size()<<endl;
			for( int j = 0; j < scafs->at( i )->size(); j++ ){
				cerr<<scafs->at( i )->at( j )->GetName()<<endl;
			}
			cerr<<endl;
			*/

			// FIXME: before calculating gaps, remove all the removable repeats and update the related edges
			// fixed
			//RemoveRepeats( scafs->at( i ) );
			CalculateGap( scafs->at( i ), scafs->at( i )->size() - 1 );
		}

		// set left and right position for each contig
		double sumLength = 0;
		double sumLengthWithGap = 0;
		for( int j = 0; j < (int) scafs->at( i )->size(); j++ ){
			scafs->at( i )->at( j )->SetRightDis( sumLength );
			scafs->at( i )->at( j )->SetRightDisWithGap( sumLengthWithGap );
			sumLength += scafs->at( i )->at( j )->GetLength() - scafs->at( i )->at( j )->GetRemovedDis();
			sumLengthWithGap += scafs->at( i ) ->at( j )->GetLength()  - scafs->at( i )->at( j )->GetRemovedDis() 
				+ scafs->at( i ) ->at( j )->GetGap();
		}
		sumLength = 0;
		sumLengthWithGap = 0;
		for( int j = scafs->at( i )->size() - 1; j >= 1; j-- ){
			scafs->at( i )->at( j )->SetLeftDis( sumLength );
			scafs->at( i )->at( j )->SetLeftDisWithGap( sumLengthWithGap );
			sumLength += scafs->at( i )->at( j )->GetLength() - scafs->at( i )->at( j )->GetRemovedDis();
			sumLengthWithGap += scafs->at( i )->at( j )->GetLength() - scafs->at( i )->at( j )->GetRemovedDis() 
				+ scafs->at( i )->at( j - 1 )->GetGap();
		}
		scafs->at( i )->at( 0 )->SetLeftDis( sumLength );
		scafs->at( i )->at( 0 )->SetLeftDisWithGap( sumLengthWithGap );
	}
	//cout<<"finish calculating gap\n";

	//fprintf( logFile, "finish calculating gaps\n" );
	//fflush( logFile );
	
	// generate scaffold
	for( int i = 0; i < (int) scafs->size(); i++ ){
		//cout<<"Handle "<<i<<" scaf\n";
		//fprintf( logFile, "Handle %d scaffold\n", i );
		//fflush( logFile );

		vector<Contig*> *currentScaf = scafs->at( i );
		int scaffoldID = currentScaf->at( 0 )->GetScafID();
		string scaffoldResult = "";
		double scaffoldLength = 0;
		double scaffoldLengthWithGap = 0;
		double scaffoldCov = 0;		// the average coverage of scaffolds
		//double endScaffoldPos = currentScaf->at( 0 )->GetStartPosition() + currentScaf->at( 0 )->GetLength() - currentScaf->at( 0 )->GetRemovedDis()- 1;
		//double endScaffoldPosWithGap = currentScaf->at( 0 )->GetStartPositionWithGap() + currentScaf->at( 0 )->GetLengthWithGap() - currentScaf->at( 0 )->GetRemovedDis() - 1;
		for( int j = 0; j < (int) currentScaf->size(); j++ ){
			//cout<<"\tgenerating "<<j<<" string\n";
			// generate the scaffold string
			Contig *addedContig = currentScaf->at( j );
			if( addedContig->NeedRemove() ){
				// if this repeat can be removed, ignore it
				continue;
			}

			/*if( j != 0 )
				scaffoldResult = "\n" + scaffoldResult;
			*/

			//cout<<"\tget contig correctly\n";
			//cout<<"\t"<<addedContig->GetName()<<endl;

			// FIXME: remove the adjacent same repeat
			//scaffoldResult = addedContig->GenScaffoldString( addedContig->GetOri() ) 
			//	+ itos( addedContig->GetGap() ) + scaffoldResult;
			//cerr<<"add string to scaffold\n";
			scaffoldResult = AddScaffoldString( addedContig->GenScaffoldString( addedContig->GetOri() ) + itos( addedContig->GetGap() ), scaffoldResult );
			scaffoldLength += addedContig->GetLength() - addedContig->GetRemovedDis();
			scaffoldLengthWithGap += addedContig->GetLengthWithGap() - addedContig->GetRemovedDis() + addedContig->GetGap();
			scaffoldCov += addedContig->GetCov() * (addedContig->GetLength() - addedContig->GetRemovedDis() );

			// calculate the left and right distance of contig
			//double leftDis = addedContig->GetStartPosition() - 1;
			//double rightDis = endScaffoldPos - addedContig->GetStartPosition() - addedContig->GetLength() - addedContig->GetRemovedDis() + 1;
			//double leftDisWithGap = addedContig->GetStartPositionWithGap() - 1;
			//double rightDisWithGap = endScaffoldPosWithGap - addedContig->GetStartPositionWithGap() - addedContig->GetLengthWithGap() 
				//+ addedContig->GetRemovedDis() + 1;

			//addedContig->SetLeftDis( leftDis );
			//addedContig->SetRightDis( rightDis );
			//addedContig->SetLeftDisWithGap( leftDisWithGap );
			//addedContig->SetRightDisWithGap( rightDisWithGap );
		}

		results[ scaffoldID ] = new ScaffoldResult();
		results[ scaffoldID ]->SetScaffoldString( scaffoldResult );
		results[ scaffoldID ]->SetLength( scaffoldLength );
		results[ scaffoldID ]->SetLengthWithGap( scaffoldLengthWithGap );
		scaffoldCov = scaffoldCov / scaffoldLength;
		results[ scaffoldID ]->SetCov( scaffoldCov );

		//cerr<<"scaffolds with gap: \n";
		//cerr<<scaffoldResult<<endl<<endl;
		//cout<<"handle this scaffold correctly\n";
		//fprintf( logFile, "handle this scaffold correctly\n" );
		//fflush( logFile );
	}
	//fprintf( logFile, "finish generating results\n" );
	//fflush( logFile );
	
}

// add scaffold string
// addtach the new string before old string
string opera::AddScaffoldString( string newString, string oldString ){
	if( oldString == "" ){
		//cerr<<"only new string: "<<endl<<newString<<endl;
		return newString;
	}

	vector<string> tempNew;
	vector<string> tempOld;
	
	Split( newString, "\n", &tempNew );
	Split( oldString, "\n", &tempOld );

	vector<string> newEnd;
	vector<string> oldStart;
	Split( tempNew.at( tempNew.size() - 1 ), "\t", &newEnd );
	Split( tempOld.at( 0 ), "\t", &oldStart );
	
	//cerr<<"new String: "<<endl<<newString<<endl;
	//cerr<<"old string: "<<endl<<oldString<<endl;
	//cerr<<"end of new string: "<<newEnd.at( 0 )<<"\t"<<newEnd.at( 1 )<<endl;
	//cerr<<"start of old string: "<<oldStart.at( 0 )<<"\t"<<oldStart.at( 1 )<<endl;
	
	string result = newString;
	if( newEnd.at( 0 ) == oldStart.at( 0 ) && newEnd.at( 1 ) == oldStart.at( 1 ) ){
		//cerr<<"find same repeat"<<newEnd.at( 0 )<<"\n";
		for( int i = 1; i < (int) tempOld.size(); i++ ){
			result += "\n" + tempOld.at( i );
		}
	}
	else{
		result += "\n" + oldString;
	}
	
	//cerr<<"combined string: "<<endl<<result<<endl;
	return result;
}

// calculate the gap sizes
void opera::CalculateGap( vector<Contig*> *contigs, int numberOfGaps ){
	/*cerr<<"calculate gaps\n";
	for( int i = 0; i < contigs->size(); i++ ){
		cerr<<contigs->at( i )->GetName()<<endl;
	}
	cerr<<endl;*/


	// initialize gap sizes
	Matrix<double> G, CE, CI;
	Vector<double> g0, ce0, ci0, x;
	//int n, m, p;
	//double sum = 0.0;

	G.resize( numberOfGaps, numberOfGaps );
	g0.resize( numberOfGaps );
	for (int i = 0; i < numberOfGaps; i++){
		for (int j = 0; j < numberOfGaps; j++)
			G[i][j] = 0;
		g0[ i ] = 0;
	}

	CE.resize( numberOfGaps, 0 );
	CI.resize( numberOfGaps, 0 );
	ce0.resize( 0 );
	ci0.resize( 0 );
	x.resize( numberOfGaps );

	// traverse each contig
#ifdef LOG
	//fprintf( logFile, "\nnew scaffolds\n" );
#endif

	string constraintsStr = ""; 
	int numberOfEdges = 0;      // the number of all concordant edges
	set<PET*> *visitedEdges = new set<PET*>;
	for( vector<Contig*>::iterator contigIter = contigs->begin(); contigIter!= contigs->end(); contigIter++ ){
		Contig *contig = *contigIter;
		//cerr<<"calculate gap: "<<contig->GetName()<<endl;
#ifdef LOG
		//fprintf( logFile, "%s\t%d\t%.0f\n", contig->GetName().c_str(), contig->GetOri(), contig->GetLength() );
#endif
		// update parameters
		list<PET*> *edges = new list<PET*>;
		// FIX ME: original edges do not have the pointer to the new repeat node
		//if( !Configure::SIMULTANEOUSLY_HANDLE_LIBRARIES ){
			edges->insert( edges->begin(), contig->GetLeftEdges()->begin(), contig->GetLeftEdges()->end() );
			edges->insert( edges->begin(), contig->GetRightEdges()->begin(), contig->GetRightEdges()->end() );
			/*}
		else{
			// this edge is super cluster, need to use original edges to calculate gap sizes
			// FIX ME: original edges do not have the pointer to the new repeat node
			list<PET*> *tempEdges = new list<PET*>;
			tempEdges->insert( tempEdges->begin(), contig->GetLeftEdges()->begin(), contig->GetLeftEdges()->end() );
			tempEdges->insert( tempEdges->begin(), contig->GetRightEdges()->begin(), contig->GetRightEdges()->end() );

			for( list<PET*>::iterator superPetIter = tempEdges->begin(); superPetIter != tempEdges->end(); superPetIter++ ){
				for( list<PET*>::iterator singlePetIter = (*superPetIter)->GetClustersOfSuperCluster()->begin();
				     singlePetIter != (*superPetIter)->GetClustersOfSuperCluster()->end(); singlePetIter++ ){
					if( (*superPetIter)->IfInSubgraph() )
						(*singlePetIter)->SetInSubgraph();
					else
						(*singlePetIter)->SetNotInSubgraph();

					if( (*superPetIter)->IfUnhappy() )
						(*singlePetIter)->SetUnhappy();
					else
						(*singlePetIter)->SetHappy();

					(*singlePetIter)->ReplaceContig( START, (*superPetIter)->GetStartContig() );
					(*singlePetIter)->ReplaceContig( END, (*superPetIter)->GetEndContig() );
					
					edges->push_back( *singlePetIter );
				}
			}

			tempEdges->clear();
			delete tempEdges;
		}
			*/

		for( list<PET*>::iterator edgeIter = edges->begin(); edgeIter != edges->end(); edgeIter++ ){
			PET *edge = *edgeIter;
			if( !edge->isInSubgraph() || edge->IfUnhappy() )
				continue;

			// check if visited this edge before
			if( visitedEdges->find( edge ) == visitedEdges->end() )
				visitedEdges->insert( edge );
			else
				continue;

			/*if( edge->GetDis() > 30000 ){
				cerr<<edge->GetOriString()<<endl;
				cerr<<"new distance: "<<edge->GetDis()<<endl;
			}
			*/

			//cerr<<edge->GetOriString()<<endl;
			//cerr<<"new distance: "<<edge->GetDis()<<endl;
			//cerr<<contig->GetName()<<endl;
			//cerr<<edge<<endl;
			//cerr<<"contigs: "<<edge->GetStartContig()<<"\t"<<edge->GetEndContig()<<endl;
			if( !contig->IsRepeat() ){
				// if the edge of a unique node connects with a repeat, get the corresponding repetitive edge
				if( edge->GetStartContig()->IsRepeat() || edge->GetEndContig()->IsRepeat() ){
					//cerr<<"it is a repetitive edge, update\n";
					if( edge->GetRepetitiveEdge() != NULL ){
						edge = edge->GetRepetitiveEdge();
						//cerr<<"get repetitive edge\n";
						//cerr<<edge->GetStartContig()<<endl;
					}
					else{
						edge = edge->GetUniqueEdge();
						//cerr<<"get unique edge\n";
					}
					//if( edge == NULL )
					//	cerr<<"repetitive edge is null\n";
				}
			}
			//cerr<<"delete edge: "<<edge->GetOriString()<<endl;

			numberOfEdges++;

			int firstPos = edge->GetStartContig()->GetScaffoldPos();
			int secondPos = edge->GetEndContig()->GetScaffoldPos();
			//cout<<edge->GetStartContig()->GetName()<<"\t"<<edge->GetStartContig()->GetScaffoldPos()<<"\t"<<edge->GetEndContig()->GetName()
			//    <<"\t"<<edge->GetEndContig()->GetScaffoldPos()<<endl;
			//cerr<<"Pos: "<<firstPos<<"\t"<<secondPos<<endl;
			if( firstPos > secondPos ){
				int temp = secondPos;
				secondPos = firstPos;
				firstPos = temp;
			}

			double sumLength = 0;		// the total contig length in the middle
			for( int id = firstPos + 1; id < secondPos; id++ ){
				sumLength += contigs->at( id )->GetLength();
				//sumLength += contigs->at( id )->GetLengthWithGap();
			}
			constraintsStr += itos( firstPos ) + "\t" + itos( secondPos ) + "\t" 
				+ itos( edge->GetDisWithGap() + Configure::STD_TIMES * edge->GetStd() - sumLength ) + "\n";
			//constraintsStr += itos( firstPos ) + "\t" + itos( secondPos ) + "\t" 
			//	+ itos( edge->GetDis() + Configure::STD_TIMES * edge->GetStd() - sumLength ) + "\n";

			// update related gap coefficients
			for( int gapID = firstPos; gapID < secondPos; gapID++ ){
				g0[ gapID ] += 2.0 * ( sumLength - (double) edge->GetDisWithGap() ) 
					/ ( (double)edge->GetStd() * (double)edge->GetStd() );
				//g0[ gapID ] += 2.0 * ( sumLength - (double) edge->GetDis() ) 
				//	/ ( (double)edge->GetStd() * (double)edge->GetStd() );
				G[ gapID ][ gapID ] += 2.0 / (double)( edge->GetStd() * edge->GetStd() );
				for( int gapID2 = firstPos; gapID2 < secondPos; gapID2++ )
				{
					if( gapID != gapID2 )
						G[ gapID ][ gapID2 ] += 2.0 / (double)( edge->GetStd() * edge->GetStd() );
				}
			}
			//cerr<<"update gap coefficient ok\n";
		}

		delete edges;
	}
	visitedEdges->clear();
	delete visitedEdges;

	/*cerr<<"G:\n";
	for( int i = 0; i < numberOfGaps; i++ ){
		for( int j = 0; j < numberOfGaps; j++ )
			cerr<<G[ i ][ j ]<<"\t";
		cerr<<endl;
	}
	cerr<<endl;*/

	// add constraints: gi>=-kmer
	CI.resize( numberOfGaps, numberOfGaps + numberOfEdges );
	ci0.resize( numberOfGaps + numberOfEdges );
	for( int i = 0; i < numberOfGaps; i++ ){
		ci0[ i ] = Configure::KMER;
		for( int j = 0; j < numberOfGaps; j++ )
		{
			if( i == j )
				CI[ i ][ j ] = 1;
			else
				CI[ i ][ j ] = 0;
		}
	}

	vector<string>* cons = new vector<string>;
	Split( constraintsStr, "\n", cons );
	for( int i = numberOfGaps; i < numberOfGaps + numberOfEdges; i++ ){
		vector<string>* line = new vector<string>;
		Split( cons->at( i - numberOfGaps ), "\t", line );
		int start = atoi( line->at( 0 ).c_str() );
		int end = atoi( line->at( 1 ).c_str() );
		int value = atoi( line->at( 2 ).c_str() );
		ci0[ i ] = value;
		for( int j = 0; j < numberOfGaps; j++ )
		{
			if( j >= start && j < end )
				CI[ j ][ i ] = -1;
			else
				CI[ j ][ i ] = 0;
		}
		line->clear();
		delete line;
	}
	cons->clear();
	delete cons;
		
	solve_quadprog(G, g0, CE, ce0, CI, ci0, x);

	// set the gap sizes to contigs
	contigs->at( 0 )->SetGap( 0 );
	for( int id = 1; id < (int) contigs->size(); id++ ){
		//if( x[ id - 1 ] > 0 )
			contigs->at( id )->SetGap( (int)x[ id - 1 ] );
		//else
		//	contigs->at( id )->SetGap( 1 );
	}

	// calculate position of each contig with gap
	double pos = 1;
        contigs->at( contigs->size() - 1 )->SetStartPositionWithGap( 1 );
	pos += contigs->at( contigs->size() - 1 )->GetLengthWithGap() + contigs->at( contigs->size() - 1 )->GetGap();
	for( int id = contigs->size() - 2; id >= 0; id-- ){
		contigs->at( id )->SetStartPositionWithGap( pos );
		//cout<<contigs->at( id )->GetStartPosition()<<"====="<<contigs->at( id )->GetStartPositionWithGap()<<endl;
		pos += contigs->at( id )->GetLengthWithGap() + contigs->at( id )->GetGap();
	}
}

// output scaffolds
// return -1 if fails
int opera::OutputScaffold( string fileName ){
	return m_graph->OutputScaffolds( fileName );
}

// output unhappy edges
int opera::OutputUnhappyEdges( string fileName ){
	ofstream unhappyWriter( fileName.c_str() );
	if( unhappyWriter == NULL ){
		cout<<"ERROR: Cannot open "<<fileName<<" file"<<endl;
		return -1;
	}

	string head = "First Contig\tOrientation\tSecond Contig\tOrientation\tDistance\tStandard Deviation\tSize\n";
	unhappyWriter.write( head.c_str(), head.length() );
	for( list<string>::iterator iter = m_unhappyEdgeString->begin(); iter != m_unhappyEdgeString->end();
		iter++ ){
			string line = (*iter) + "\n";
			unhappyWriter.write( line.c_str(), line.length() );
	}

	unhappyWriter.close();
	return 1;
}

// print scaffold result
void opera::PrintScaffoldResult( ScaffoldResult **results, int number ){
	for( int i = 0; i < number; i++ ){
		cout<<"Scaffold "<<i<<endl;
		cout<<results[ i ]->GetScaffoldString()<<endl;
	}
}

// Get statistics of assembly
string opera::GetStats(){
	return m_stats;
}

// sort all scaffold according to length
// return -1 if failed
int opera::SortScaffold(){
	multiset<FinalScaffold*, great_length> *scaffolds = new multiset<FinalScaffold*, great_length>;
	double totalLength = 0;

	// read small contigs
	//ReadContigFile( Configure::OUTPUT_FOLDER + Configure::DIRECTORY_SIGN + "smallContigs", scaffolds, totalLength );
	//cerr<<"After read small contig, total Length is "<<totalLength<<endl;
	
	// read repeats
	//ReadContigFile( Configure::OUTPUT_FOLDER + Configure::DIRECTORY_SIGN + "repeatContigs", scaffolds, totalLength );
	//cerr<<"After read repeat: "<<totalLength<<endl;

	// read scaffold
	int nonSingleton = 0;
	ReadScaffoldFile( Configure::OUTPUT_FOLDER + Configure::DIRECTORY_SIGN + "scaffolds.scaf", scaffolds, totalLength, nonSingleton );
	//cerr<<"After read scaffold: "<<totalLength<<endl; 

	// calcualte N50
	double N50 = 0;
	multiset<FinalScaffold*, great_length>::iterator iter = scaffolds->begin();
	while( iter != scaffolds->end() ){
		N50 += (*iter)->GetLength();
		if( N50 > 0.5 * totalLength ){
			N50 = (*iter)->GetLength();
			break;
		}
		iter++;
	}
	iter = scaffolds->begin();
	string fileName = Configure::OUTPUT_FOLDER + Configure::DIRECTORY_SIGN + "statistics";

	// write in another way
	m_stats = "\nScaffold statistics:\n";
	m_stats += "\tN50: ";
	//cout<<endl<<"Scaffold statistics:"<<endl;
	//cout<<"\tN50: ";
	m_stats += PrintNumber( N50 );
	m_stats += "\n";
	//cout<<"\tTotal length: ";
	m_stats += "\tTotal length: ";
	m_stats += PrintNumber( totalLength );
	m_stats += "\n";
	//cout<<endl;
	//cout<<"\tLongest scaffold: ";
	m_stats += "\tLongest scaffold: ";
	m_stats += PrintNumber( (*iter)->GetLength() );
	m_stats += "\n";
	//cout<<endl;

	FILE *contigWriter = fopen( fileName.c_str(), "w" );
	fprintf( contigWriter, "N50: " );
	PrintNumberToFile( contigWriter, N50 );
	fprintf( contigWriter, "\n" );
	fprintf( contigWriter, "Total length of assembly: " );
	PrintNumberToFile( contigWriter, totalLength );
	fprintf( contigWriter, "\n" );
	fprintf( contigWriter, "Number of non-singleton scaffolds: %'i \n", nonSingleton );
	fprintf( contigWriter, "Length of the longest scaffold: " );
	PrintNumberToFile( contigWriter, (*iter)->GetLength() );
	fprintf( contigWriter, "\n" );
	fclose( contigWriter );

	// read original fasta file
	map<string, string> *contigs = new map<string, string>;
	ReadContigFasta( Configure::CONTIG_FILE, contigs );

	// output results
	fileName = Configure::OUTPUT_FOLDER + Configure::DIRECTORY_SIGN + "scaffoldSeq.fasta";
	ofstream resultWriter( fileName.c_str() );

	if( resultWriter == NULL ){
		cout<<"ERROR: Cannot open "<<fileName<<" file"<<endl;
		return -1;
	}
	iter = scaffolds->begin();
	vector<string> *scafVec = new vector<string>;
	vector<string> *scafLine = new vector<string>;
	while( iter != scaffolds->end() ){
		string resultString = ">" + (*iter)->GetName();
		resultWriter.write( resultString.c_str(), resultString.length() );

		if( !(*iter)->IsScaffold() ){
			// small contigs and repeats
			char tempName[500];
			sprintf( tempName, "\tlength: %.0f\tcov: %.1f\n", (*iter)->GetLength(), (*iter)->GetCov() );
			string seq( tempName );
			resultWriter.write( seq.c_str(), seq.length() );

			seq = (*contigs)[ (*iter)->GetName() ];
			resultWriter.write( seq.c_str(), seq.length() );
		}
		else{
			// scaffolds
			resultString = "\n";
			resultWriter.write( resultString.c_str(), resultString.length() );

			string scafString = (*iter)->GetScaffold();
			Split( scafString, "\n", scafVec );
			string content = "";
			for( int i = 0; i < (int) scafVec->size(); i++ ){
				scafLine->clear();
				Split( (*scafVec)[ i ], "\t", scafLine );
				content = (*contigs)[ (*scafLine)[ 0 ] ];
				// output each contig
				if( (*scafLine)[ 1 ] == "BE" ){
					// output "+"
					resultWriter.write( content.c_str(), (*contigs)[ (*scafLine)[ 0 ] ].length() );
				}
				else{
					// output "-"
					string reverseSeq = "";
					for ( string::reverse_iterator reverseIter = content.rbegin(); reverseIter != content.rend(); reverseIter++ )
					{
						if( *reverseIter == 'A' || *reverseIter == 'a')
							reverseSeq.append( "T" );
						else if( *reverseIter == 'C' || *reverseIter == 'c')
							reverseSeq.append( "G" );
						else if( *reverseIter == 'G' || *reverseIter == 'g' )
							reverseSeq.append( "C" );
						else if( *reverseIter == 'T' || *reverseIter == 't' )
							reverseSeq.append( "A" );
						else
							reverseSeq.append( "N" );
					}
					resultWriter.write( reverseSeq.c_str(), reverseSeq.size() );
				}
				int gapSize = atoi( (*scafLine)[ 3 ].c_str() );
				// add the gaps except for the last contig
				if( gapSize < 10 && i != (int) scafVec->size() -1 )
					gapSize = 10;
				string gap = "";
				for( int num = 0; num < gapSize; num++ )
					gap.append( "N" );
				resultWriter.write( gap.c_str(), gap.size() );
			}
			scafVec->clear();
		}
		string newLine = "\n";
		resultWriter.write( newLine.c_str(), newLine.length() );
		iter++;
	}

	resultWriter.close();

	// clear
	iter = scaffolds->begin();
	while( iter != scaffolds->end() ){
		//cout<<(*iter)->GetName()<<"\t"<<(*iter)->GetLength()<<endl;
		delete *iter;
		scaffolds->erase( iter++ );
	}
	delete scaffolds;

	contigs->clear();
	delete contigs;
	scafLine->clear();
	delete scafLine;
	scafVec->clear();
	delete scafVec;
	return 0;
}

// read result contig files
int opera::ReadContigFile( string fileName, multiset<FinalScaffold*, great_length> *scaffolds, double &totalLength )
{
	ifstream contigReader( fileName.c_str() );

	if( contigReader == NULL )
	{
		cout<<"ERROR: Cannot open "<<fileName<<" file"<<endl;
		return -1;
	}

	string line;
	vector<string> *contents = new vector<string>;
	getline( contigReader, line );		// remove the first title line
	while( getline( contigReader, line ) ){
		Split( line, "\t", contents );
		FinalScaffold *newFS = new FinalScaffold( (*contents)[ 0 ], atoi( (*contents)[ 1 ].c_str() ), "", atoi( (*contents)[ 2 ].c_str() ) );
		totalLength += newFS->GetLength();
		scaffolds->insert( newFS );
	}

	delete contents;
	contigReader.close();
	return 0;
}

// read scaffolds file
int opera::ReadScaffoldFile( string fileName, multiset<FinalScaffold*, great_length> *scaffolds, double &totalLength, int &nonSingleton )
{
	ifstream contigReader( fileName.c_str() );

	if( contigReader == NULL )
	{
		cout<<"ERROR: Cannot open "<<fileName<<" file"<<endl;
		return -1;
	}

	string line;
	vector<string> *contents = new vector<string>;
	bool isFirst = true;
	int length = 0;
	string name = "";
	string scaffold = "";
	int numOfContig = 0;
	while( getline( contigReader, line ) ){
		if( line.at( 0 ) == '>' ){
			// save previous scaffold
			if( isFirst )
				isFirst = false;
			else{
				FinalScaffold *newFS = new FinalScaffold( name, length, scaffold );
				totalLength += length;
				scaffolds->insert( newFS );
				if( numOfContig > 1 )
					nonSingleton++;
			}

			// start a new scaffold
			length = 0;
			name = line.substr( 1 );
			scaffold = "";
			numOfContig = 0;
		}
		else{
			Split( line, "\t", contents );
			scaffold += line + "\n";
			length += atoi( (*contents)[ 2 ].c_str() );
			numOfContig++;
		}
	}
	// save the last scaffold
	FinalScaffold *newFS = new FinalScaffold( name, length, scaffold );
	totalLength += length;
	scaffolds->insert( newFS );
	if( numOfContig > 1 )
		nonSingleton++;

	contigReader.close();
	delete contents;
	return 0;
}

// read contig fasta file
int opera::ReadContigFasta( string fileName, map<string, string> *contigs )
{
	ifstream contigReader( fileName.c_str() );

	if( contigReader == NULL )
	{
		cout<<"ERROR: Cannot open "<<fileName<<" file"<<endl;
		return -1;
	}

	string line;
	string seq = "";
	string name = "";
	vector<string> *column = new vector<string>;
	while( getline( contigReader, line ) ){
		if( line.length() == 0 )
			continue;
		if( line.at( 0 ) == '>' ){
			// save previous contig
			if( seq.length() > 0 ){
				contigs->insert( pair<string, string>( name, seq ) );
			}

			// start new contig
			Split( line, " ", column );
			Split( (*column)[ 0 ], "\t", column );
			name = (*column)[ 0 ].substr( 1, (*column)[ 0 ].length() - 1 );
			seq = "";
		}
		else{
			seq.append( line );
		}
	}
	column->clear();
	delete column;
	contigReader.close();
	// save last contig
	contigs->insert( pair<string, string>( name, seq ) );
	return 0;
}

// check contig file format: velvet, soap or normal fasta file
int opera::CheckContigFormat()
{
	ifstream contigReader( Configure::CONTIG_FILE.c_str() );

	if( contigReader == NULL )
	{
		cout<<"ERROR: Cannot open "<<Configure::CONTIG_FILE<<" file"<<endl;
		return -1;
	}

	string line;
	getline( contigReader, line );

	if( line.find( "NODE_" ) != string::npos && line.find( "_length_" ) != string::npos && line.find( "_cov_" ) != string::npos ){
	//if( line.find( "contig_" ) != string::npos && line.find( "_length_" ) != string::npos && line.find( "_cov_" ) != string::npos ){
		// velvet format
		Configure::FILE_TYPE = VELVET;
		//Configure::FILTER_REPEAT = true;
	}
	else{
		vector<string> *column = new vector<string>;
		Split( line, " ", column );
		if( ( column->size() == 2 ) && IsNumber( (*column)[ 1 ] ) ){  //">scaffold" ){
			// soap scaffold format
			//cout<<"soap format\n";
			Configure::FILE_TYPE = SOAP;
			//Configure::FILTER_REPEAT = true;
		}
		else{
			if( line.find( "_tip_" ) != string::npos && line.find( " length " ) != string::npos && line.find( " cvg_" ) != string::npos ){
				// soap contig format
				//cout<<"soap contig format\n";
				Configure::FILE_TYPE = SOAP_CONTIG;
				//Configure::FILTER_REPEAT = true;
			}
			else if( line.find( "length:" ) != string::npos && line.find( "cov:" ) != string::npos ){
				// opera scaf format
				Configure::FILE_TYPE = OPERA_SCAF;
			}
			else{
				// normal fasta format
				//cout<<"normal\n";
				Configure::FILE_TYPE = NORMAL;
				//Configure::FILTER_REPEAT = true;
			}
		}
		delete column;
	}
	contigReader.close();

	return 1;
}

void opera::OpenLogFile()
{
	logFile = fopen( ( Configure::OUTPUT_FOLDER + "log" ).c_str(), "w" );
	m_discardEdgeFile = fopen( (Configure::OUTPUT_FOLDER + "discardEdges").c_str(), "w" );
}

// move the edges in first library into graph
void opera::MoveEdgesToGraph( PetLibrary *lib ){
	m_graph->ClearEdge();

	// add edge to graph
	list<PET*> *edges = lib->GetEdges();
	list<PET*>::iterator edgeIter = edges->begin();
	//cout<<edges->size()<<" edges in this library\n";
	//cerr<<lib->GetFileName()<<endl;
	int tempNum = 0;
	while( edgeIter != edges->end() ){
		PET *edge = *edgeIter;

		//if( edge->GetStartContig()->GetName() == "10639727" || edge->GetEndContig()->GetName() == "10639727" ){
		//	cerr<<edge->ToString()<<endl;
		//}


		bool isBigContig = true;
		if( lib->GetFileName() != "super library (small libraries)" ){
			//cerr<<"not small library edge\n";
			if( edge->GetStartContig()->GetLengthWithGap() < Configure::CONTIG_SIZE_THERSHOLD 
			    || edge->GetEndContig()->GetLengthWithGap() < Configure::CONTIG_SIZE_THERSHOLD ){
				isBigContig = false;
				//cerr<<edge->GetOriString()<<endl;
			}

			//if( edge->GetStartContig()->GetName() == "10639727" || edge->GetEndContig()->GetName() == "10639727" ){
			//	cerr<<edge->ToString()<<endl;
			//}
		}

		/*if( edge->GetOriginalStartContig()->GetName() == "1601" && edge->GetOriginalEndContig()->GetName() == "1613" ){
			cerr<<"find edge: "<<edge->GetOriString()<<"\tDis: "<<edge->GetDis()<<endl;
		}
		*/
		// FIXME: check edges to see if two repeats are adjacent to each other
		// FIXED
		m_graph->CheckEdgesForAdjacentRepeat( edge );

		/*
		if( edge->GetOriginalStartContig()->GetName() == "1601" && edge->GetOriginalEndContig()->GetName() == "1613" ){
			cerr<<"find edge(after modification): "<<edge->GetOriString()<<"\tDis: "<<edge->GetDis()<<endl;
		}
		*/

		bool ifDeleteEdge = false;
		if( isBigContig && edge->GetDis() + Configure::STD_TIMES * edge->GetStd() >= -Configure::KMER ){
			// a valid edge, save to graph
			m_graph->AddEdge( edge );
			tempNum++;
			//cout<<"New edge\t"<<edge->GetOriString()<<endl;
		}
		else
		{
			m_unhappyEdgeString->push_back( edge->GetOriString() );
			ifDeleteEdge = true;
			//cout<<"Invalid distance: "<<edge->GetDis()<<":\t"<<edge->GetOriString()<<endl;
		}

		// remove edge from unused edges list in contigs
		edge->GetStartContig()->RemoveEdgeMultiLib( edge, START );
		edge->GetEndContig()->RemoveEdgeMultiLib( edge, END );

		/*		if( edge->GetStartContig()->GetName() == "1241" && edge->GetEndContig()->GetName() == "1455" ){
			cerr<<"FIND EDGE: "<<edge->GetOriString()<<endl<<edge<<endl;
		}
		*/

		edgeIter = edges->erase( edgeIter );

		if( ifDeleteEdge )
			delete edge;
	}

	// find 1603 contig
	//m_graph->CheckNumberOfEdges( "1607" );
}

// analyze assembly file to get kmer size
bool opera::GetKmer( string contigFileName, int type ){
	// get the directory path of assembly
	vector<string> *names = new vector<string>;
	string path = "";
	Split( contigFileName, "/", names );
	if( names->size() != 1 ){
		path = contigFileName.substr( 0, contigFileName.length() - names->at( names->size() - 1 ).length() );
	}
	string fileName = names->at( names->size() - 1 );
	names->clear();
	delete names;

	// if the assembly was produced by velvet, try to find LastGraph file
	if( type == VELVET ){
		string graphFile = path + "LastGraph";
		ifstream graphReader( graphFile.c_str() );

		if( graphReader == NULL ){
			return false;
		}

		string line;
		getline( graphReader, line );
		vector<string> *temp = new vector<string>;
		Split( line, "\t", temp );
		Configure::KMER = atoi( temp->at( 2 ).c_str() );
		temp->clear();
		delete temp;

		graphReader.close();
		return true;
	}
		
	// if the assembly was produced by SOAPdenovo, try to find preGraphBasic file
	if( type == SOAP || type == SOAP_CONTIG ){
		// get the prefix
		vector<string> *temp = new vector<string>;
		Split( fileName, ".", temp );
		string graphFile = path + temp->at( 0 );
		graphFile += ".preGraphBasic";

		ifstream graphReader( graphFile.c_str() );

		if( graphReader == NULL ){
			temp->clear();
			delete temp;
			return false;
		}

		string line;
		getline( graphReader, line );
		temp->clear();
		Split( line, " ", temp );
		Configure::KMER = atoi( temp->at( 3 ).c_str() );
		temp->clear();
		delete temp;

		graphReader.close();
		return true;
	}

	return false;
}

// print the parameter into a file
int opera::PrintParameter( string fileName ){
	ofstream paraWriter( fileName.c_str() );
	
	if( paraWriter == NULL ){
		cout<<"ERROR: Cannot open "<<fileName<<" file\n";
		return -1;
	}

	string resultString = "";
	resultString = "kmer\t" + itos( Configure::KMER ) + "\n";
	//resultString += "threshold\t" + itos( Configure::CLUSTER_THRESHOLD ) + "\n";
	
	paraWriter.write( resultString.c_str(), resultString.length() );

	paraWriter.close();
	return 1;
}

// output conflicting edges to file
int opera::OutputConflictingEdges( string fileName ){
	ofstream confWriter( fileName.c_str() );
	
	if( confWriter == NULL ){
		cout<<"ERROR: Cannot open "<<fileName<<" file\n";
		return -1;
	}

	string resultString = "First Contig\tOrientation\tSecond Contig\tOrientation\tDistance\tStandard Deviation\tSize\n";
	for( list<string>::iterator confIter = m_conflictingEdges->begin(); confIter != m_conflictingEdges->end(); confIter++ ){
		resultString.append( *confIter );
		resultString.append( "\n" );
	}
	
	confWriter.write( resultString.c_str(), resultString.length() );

	confWriter.close();
	return 1;
}

// check if a repeat can be removed from the active region without introducing any discordant edges
// c: the repeat about to be removed
// ori: the orientation of c
// pos: the position of c in Active Region (from the tail)
// p: the current partial scaffold
/*bool opera::BringDiscordantEdgeByRemoving( Contig *c, int ori, int pos, PartialScaffold *p ){
	     // find the region surrounding c (bounded by the upperbound)
	int upperbound = Configure::LIB_MEAN + Configure::STD_TIMES * Configure::LIB_STD;
	list<SUBCONTIG*> wholeRegion;

	// find the region behind c
	list<SUBCONTIG*> regionAfterC;
	int counter = m_activeRegion->size();
	int length = 0;
	for( list<Contig*>::iterator arIter = m_activeRegion->begin(); arIter != m_activeRegion->end(); arIter++ ){
		if( counter >= pos ){
			continue;
		}
		else{
			//regionAfterC.push_back( *arIter );
			// add subcontigs
			bool foundSameContig = false;
			list<SUBCONTIG*> *subcontigs = (*arIter)->GetSubcontigs();
			for( list<SUBCONTIG*>::iterator subIter = subcontigs->begin(); subIter != subcontigs->end(); subIter++ ){
				wholeRegion.push_back( *subIter );
				regionAfterC.push_back( *subIter );
				length += (*subIter)->m_length;

				if( (*subIter)->m_contigName == c->GetName() && (*subIter)->m_ori == ori ){
					foundSameContig = true;
					break;
				}
			}

			if( length > upperbound || foundSameContig )
				break;
		}
	}

	// find the region before c
	list<SUBCONTIG*> regionBeforeC;
	counter = 1;
	length = 0;
	for( list<Contig*>::reverse_iterator arIter = m_activeRegion->rbegin(); arIter != m_activeRegion->rend(); arIter++ ){
		if( counter <= pos ){
			continue;
		}
		else{
			// add subcontigs
			bool foundSameContig = false;
			list<SUBCONTIG*> *subcontigs = (*arIter)->GetSubcontigs();
			
			for( list<SUBCONTIG*>::reverse_iterator subIter = subcontigs->rbegin(); subIter != subcontigs->rend(); subIter++ ){
				wholeRegion.push_front( *subIter );
				regionBeforeC.push_front( *subIter );
				length += (*subIter)->m_length;

				if( (*subIter)->m_contigName == c->GetName() && (*subIter)->m_ori == ori ){
					foundSameContig = true;
					break;
				}
			}

			//regionBeforeC.push_front( *arIter );
			if( length > upperbound || foundSameContig )
				break;
		}
	}
	
	while( length < upperbound ){
		// continue to add contigs by tracing back the scaffold
		PartialScaffold *currentP = p->GetParent();
		Contig *currentContig = currentP->GetAddedContigInAR();
		if( currentContig != NULL ){
			//regionBeforeC.push_front( currentContig );
			// add subcontigs
			bool foundSameContig = false;
			list<SUBCONTIG*> *subcontigs = currentContig->GetSubcontigs();
			for( list<SUBCONTIG*>::reverse_iterator subIter = subcontigs->rbegin(); subIter != subcontigs->rend(); subIter++ ){
				wholeRegion.push_front( *subIter );
				regionBeforeC.push_front( *subIter );
				length += (*subIter)->m_length;

				if( (*subIter)->m_contigName == c->GetName() && (*subIter)->m_ori == ori ){
					foundSameContig = true;
					break;
				}
			}

			if( length > upperbound || foundSameContig )
				break;
		}
		else{
			// reach to the beginning of the scaffold
			break;
		}
	}

	// find the concordant edges connecting with c
	list<PET*> leftEdges;      // all the edges of c connecting the nodes on the left of c
	list<PET*> rightEdges;     // all the edges of c connecting the nodes on the right of c

	list<PET*> *possibleLeftEdges, *possibleRightEdges;
	if( ori == PLUS ){
		possibleLeftEdges = c->GetLeftEdges();
		possibleRightEdges = c->GetRightEdges();
	}
	else{
		possibleLeftEdges = c->GetRightEdges();
		possibleRightEdges = c->GetLeftEdges();
	}

	// collect left edges
	CollectRelatedEdges( &leftEdges, possibleLeftEdges, c, &regionBeforeC );
	
	// collect right edges
	CollectRelatedEdges( &rightEdges, possibleRightEdges, c, &regionAfterC );

	// create new edges for this repeat
	AddRepetitiveEdges( leftEdges, rightEdges, c );

	// check if all the edges can be satisfied or not, after removing c
	list<PET*> allEdges;
	allEdges.insert( allEdges.begin(), leftEdges.begin(), leftEdges.end() );
	allEdges.insert( allEdges.begin(), rightEdges.begin(), rightEdges.end() );

	bool needToRemove = CheckEdgeSatisfaction( &allEdges, &wholeRegion, c );
	c->SetIfRemove( needToRemove );

	return true;
	
}*/

bool opera::BringDiscordantEdgeByRemoving( Contig *c, int ori, list<Contig*> erasedContigs, PartialScaffold *p, int sizeOfActiveRegion ){
	//cerr<<"check discordant edge by removing a repeat\n";
	//cerr<<"\tremoving "<<c->GetName()<<endl;

	// find the region surrounding c (bounded by the upperbound)
	int upperbound = Configure::LIB_MEAN + Configure::STD_TIMES * Configure::LIB_STD;
	//cerr<<"upperbound is "<<upperbound<<endl;
	list<SUBCONTIG*> wholeRegion;

	//cerr<<"checking repeat: "<<c->GetName()<<endl;
	// find the region behind c
	//cerr<<"find the region behind c\n";
	list<SUBCONTIG*> regionAfterC;
	int length = 0;
	for( list<Contig*>::iterator arIter = m_activeRegion->begin(); arIter != m_activeRegion->end(); arIter++ ){
		// add subcontigs
		list<SUBCONTIG*> *subcontigs = (*arIter)->GetSubcontigs();
		for( list<SUBCONTIG*>::iterator subIter = subcontigs->begin(); subIter != subcontigs->end(); subIter++ ){
			wholeRegion.push_back( *subIter );
			regionAfterC.push_back( *subIter );
			length += (*subIter)->m_length;

			if( length > upperbound + Configure::KMER)
				break;
		}

		if( length > upperbound + Configure::KMER)
			break;
	}

	/*cerr<<"number of contigs after repeat: "<<regionAfterC.size()<<endl;
	for( list<SUBCONTIG*>::iterator sIter = regionAfterC.begin(); sIter != regionAfterC.end(); sIter++ ){
		cerr<<(*sIter)->m_contigName<<"\t";
		if( (*sIter)->m_ori == PLUS )
			cerr<<"+\t";
		else
			cerr<<"-\t";
	}
	cerr<<endl;
	*/
	
	

	// find the region before c
	//cerr<<"find the region after c\n";
	list<SUBCONTIG*> regionBeforeC;
	length = 0;
	for( list<Contig*>::reverse_iterator erasedIter = erasedContigs.rbegin(); erasedIter != erasedContigs.rend(); erasedIter++ ){
		// don't add the removed repeat
		if( erasedIter == erasedContigs.rbegin() )
			continue;

		// don't consider the removed contigs
		if( (*erasedIter)->NeedRemove() )
			continue;

		// add subcontigs
		list<SUBCONTIG*> *subcontigs = (*erasedIter)->GetSubcontigs();
			
		for( list<SUBCONTIG*>::reverse_iterator subIter = subcontigs->rbegin(); subIter != subcontigs->rend(); subIter++ ){
			wholeRegion.push_front( *subIter );
			regionBeforeC.push_front( *subIter );
			length += (*subIter)->m_length;
			
			if( length > upperbound + Configure::KMER)
				break;
		}

		if( length > upperbound + Configure::KMER)
			break;
	}
	/*cerr<<"number of contigs before repeat: "<<regionBeforeC.size()<<endl;
	for( list<SUBCONTIG*>::iterator sIter = regionBeforeC.begin(); sIter != regionBeforeC.end(); sIter++ ){
		cerr<<(*sIter)->m_contigName<<"\t";
		if( (*sIter)->m_ori == PLUS )
			cerr<<"+\t";
		else
			cerr<<"-\t";
	}
	cerr<<endl;
	*/
	
	//cerr<<"continue tracking\n";
	PartialScaffold *currentP = p;
	// trace back sizeOfActiveRegion times
	for( int i = 0; i < sizeOfActiveRegion - 1; i++ ){
		currentP = currentP->GetParent();
	}

	while( length <= upperbound + Configure::KMER ){
		// continue to add contigs by tracing back the scaffold
		currentP = currentP->GetParent();
		Contig *currentContig = currentP->GetAddedContigInAR();

		if( currentContig != NULL ){
			// don't consider the removed contigs
			if( currentContig->NeedRemove() )
				continue;

			// add subcontigs
			list<SUBCONTIG*> *subcontigs = currentContig->GetSubcontigs();
			for( list<SUBCONTIG*>::reverse_iterator subIter = subcontigs->rbegin(); subIter != subcontigs->rend(); subIter++ ){
				wholeRegion.push_front( *subIter );
				regionBeforeC.push_front( *subIter );
				length += (*subIter)->m_length;

				if( length > upperbound + Configure::KMER)
					break;
			}

			if( length > upperbound + Configure::KMER)
				break;
		}
		else{
			// reach to the beginning of the scaffold
			break;
		}
	}

	/*cerr<<"number of contigs before repeat: "<<regionBeforeC.size()<<endl;
	for( list<SUBCONTIG*>::iterator sIter = regionBeforeC.begin(); sIter != regionBeforeC.end(); sIter++ ){
		cerr<<(*sIter)->m_contigName<<"\t";
		if( (*sIter)->m_ori == PLUS )
			cerr<<"+\t";
		else
			cerr<<"-\t";
	}
	cerr<<endl;
	*/
	


	// find the concordant edges connecting with c
	list<PET*> leftEdges;      // all the edges of c connecting the nodes on the left of c
	list<PET*> rightEdges;     // all the edges of c connecting the nodes on the right of c

	list<PET*> *possibleLeftEdges, *possibleRightEdges;
	if( ori == PLUS ){
		possibleLeftEdges = c->GetOriginalRepeat()->GetLeftEdges();
		possibleRightEdges = c->GetOriginalRepeat()->GetRightEdges();
	}
	else{
		possibleLeftEdges = c->GetOriginalRepeat()->GetRightEdges();
		possibleRightEdges = c->GetOriginalRepeat()->GetLeftEdges();
	}
	

	//cerr<<"possible left edges: "<<possibleLeftEdges->size()<<endl;
	//cerr<<"possible right edges: "<<possibleRightEdges->size()<<endl;

	//cerr<<"collect edges\n";
	// collect left edges
	//cerr<<"collect left edges\n";
	//CollectRelatedEdges( &leftEdges, possibleLeftEdges, c, &regionBeforeC );
	CollectRelatedLeftEdges( &leftEdges, possibleLeftEdges, c, &regionBeforeC );
	
	// collect right edges
	//cerr<<"collect right edges\n";
	//CollectRelatedEdges( &rightEdges, possibleRightEdges, c, &regionAfterC );
	CollectRelatedRightEdges( &rightEdges, possibleRightEdges, c, &regionAfterC );

	//cerr<<"left edges: "<<leftEdges.size()<<endl;
	//cerr<<"right edges: "<<rightEdges.size()<<endl;

	//cerr<<"create new edges for this repeat\n";
	// create new edges for this repeat
	AddRepetitiveEdges( leftEdges, rightEdges, c );

	// check if all the edges can be satisfied or not, after removing c
	list<PET*> allEdges;
	allEdges.insert( allEdges.begin(), leftEdges.begin(), leftEdges.end() );
	allEdges.insert( allEdges.begin(), rightEdges.begin(), rightEdges.end() );

	/*cerr<<"All related edges: \n";
	for( list<PET*>::iterator edgeIter = allEdges.begin(); edgeIter != allEdges.end(); edgeIter++ ){
		cerr<<(*edgeIter)->GetOriString()<<endl;
	}
	cerr<<endl;
	*/

	//cerr<<"check edge satisfaction\n";
	bool needToRemove = CheckEdgeSatisfaction( &allEdges, &wholeRegion, c );
	c->SetIfRemove( needToRemove );

	//cerr<<"finish check edge satisfaction\n";
	return needToRemove;
	
}

// add the repetitive edges to repeats
void opera::AddRepetitiveEdges( list<PET*> leftEdges, list<PET*> rightEdges, Contig *repeat ){
	list<PET*> allEdges;
	allEdges.insert( allEdges.end(), leftEdges.begin(), leftEdges.end() );
	allEdges.insert( allEdges.end(), rightEdges.begin(), rightEdges.end() );
	
	for( list<PET*>::iterator iter = allEdges.begin(); iter != allEdges.end(); iter++ ){
		//int repeatPos = (*iter)->GetPositionOfContig( repeat->GetOriginalRepeat() );
		int repeatOri = (*iter)->GetOrientationOfContig( repeat->GetOriginalRepeat() );
		Contig *otherContig = (*iter)->GetOtherContig( repeat->GetOriginalRepeat() );
		int otherPos = (*iter)->GetPositionOfContig( otherContig );
		int otherOri = (*iter)->GetOrientationOfContig( otherContig );

		PET *newEdge;
		//cerr<<"CREATE new repetitive edge: "<<(*iter)->GetOriString()<<endl;
		if( otherPos == START ){
			newEdge = new PET( otherContig, otherOri, repeat, repeatOri, (*iter)->GetDis(), (*iter)->GetStd(), (*iter)->GetSize(), *iter );
		}
		else{
			newEdge = new PET( repeat, repeatOri, otherContig, otherOri, (*iter)->GetDis(), (*iter)->GetStd(), (*iter)->GetSize(), *iter );
		}

		(*iter)->SetRepetitiveEdge( newEdge );
		newEdge->SetOriginalEdge( *iter );
		repeat->AddEdge( newEdge );
		(*iter)->SetUniqueEdge( NULL );
		newEdge->SetDisWithGap( (*iter)->GetDisWithGap() );
		//cerr<<(*iter)<<endl;
		//cerr<<"contigs: "<<(*iter)->GetStartContig()<<"\t"<<(*iter)->GetEndContig()<<endl;
		//cerr<<newEdge<<endl;
		//cerr<<"contigs: "<<newEdge->GetStartContig()<<"\t"<<newEdge->GetEndContig()<<endl;

		/*		if( (*iter)->GetEndContig()->GetName() == "1455" ){
			// print the edge
			for( list<PET*>::iterator tempEdgeIter = (*iter)->GetEndContig()->GetRightEdges()->begin(); 
			     tempEdgeIter != (*iter)->GetEndContig()->GetRightEdges()->end(); tempEdgeIter++ ){
				cerr<<(*tempEdgeIter)->GetOriString()<<endl;
				cerr<<(*tempEdgeIter)<<endl;
			}
		}
		*/

		//cerr<<(*iter)->GetRepetitiveEdge()->GetStartContig()->GetName()<<"\t"<<(*iter)->GetRepetitiveEdge()->GetEndContig()->GetName()<<endl;
	}
}

// add the repetitive edges to repeats
void opera::AddRepetitiveEdges( PET *edge, Contig *repeat ){
	//int repeatPos = edge->GetPositionOfContig( repeat->GetOriginalRepeat() );
	int repeatOri = edge->GetOrientationOfContig( repeat->GetOriginalRepeat() );
	Contig *otherContig = edge->GetOtherContig( repeat->GetOriginalRepeat() );
	int otherPos = edge->GetPositionOfContig( otherContig );
	int otherOri = edge->GetOrientationOfContig( otherContig );

	PET *newEdge;
	//cerr<<"CREATE new repetitive edge: "<<edge->GetOriString()<<endl;
	if( otherPos == START ){
		newEdge = new PET( otherContig, otherOri, repeat, repeatOri, edge->GetDis(), edge->GetStd(), edge->GetSize(), edge );
	}
	else{
		newEdge = new PET( repeat, repeatOri, otherContig, otherOri, edge->GetDis(), edge->GetStd(), edge->GetSize(), edge );
	}

	edge->SetRepetitiveEdge( newEdge );
	newEdge->SetOriginalEdge( edge );
	repeat->AddEdge( newEdge );
	edge->SetUniqueEdge( NULL );
	//cerr<<(*iter)<<endl;
	//cerr<<"contigs: "<<(*iter)->GetStartContig()<<"\t"<<(*iter)->GetEndContig()<<endl;
	//cerr<<newEdge<<endl;
	//cerr<<"contigs: "<<newEdge->GetStartContig()<<"\t"<<newEdge->GetEndContig()<<endl;
	
	/*		if( (*iter)->GetEndContig()->GetName() == "1455" ){
	// print the edge
	for( list<PET*>::iterator tempEdgeIter = (*iter)->GetEndContig()->GetRightEdges()->begin(); 
	tempEdgeIter != (*iter)->GetEndContig()->GetRightEdges()->end(); tempEdgeIter++ ){
	cerr<<(*tempEdgeIter)->GetOriString()<<endl;
	cerr<<(*tempEdgeIter)<<endl;
	}
	}
	*/
	
	//cerr<<(*iter)->GetRepetitiveEdge()->GetStartContig()->GetName()<<"\t"<<(*iter)->GetRepetitiveEdge()->GetEndContig()->GetName()<<endl;
}

// check if all the edges can be satisfied
// if the number of discordant edges does not change, return true; else return false
/*
bool opera::CheckEdgeSatisfaction( list<PET*> *collectedEdges, list<SUBCONTIG*> *region, Contig *repeat ){

	// FIXME: check after removing this repeat, if the number of discordant edges will increase

	// create a map of all subcontigs
	// map of forward order
	map<string, list<SUBCONTIG*>::iterator> mapOfContigs;
	// map of reverse order
	map<string, list<SUBCONTIG*>::reverse_iterator> reverseMapOfContigs;
	
	for( list<SUBCONTIG*>::iterator subIter = region->begin(); subIter != region->end(); subIter++ ){
		SUBCONTIG *currentContig = *subIter;
		if( !currentContig->m_isRepeat ){
			// add unique contigs
			mapOfContigs.insert( pair<string, list<SUBCONTIG*>::iterator> (currentContig->m_contigName, subIter) );
		}
	}

	for( list<SUBCONTIG*>::reverse_iterator subIter = region->rbegin(); subIter != region->rend(); subIter++ ){
		SUBCONTIG *currentContig = *subIter;
		if( !currentContig->m_isRepeat ){
			// add unique contigs
			reverseMapOfContigs.insert( pair<string, list<SUBCONTIG*>::reverse_iterator> (currentContig->m_contigName, subIter) );
		}
	}

	// count the number of new discordant edges
	int numOfNewDiscordantEdges = 0;
	SUBCONTIG *firstContig = *(region->begin());
	SUBCONTIG *lastContig = *(region->rbegin());
	// for each edge, check if it is can be satisfied
	for( list<PET*>::iterator edgeIter = collectedEdges->begin(); edgeIter != collectedEdges->end(); edgeIter++ ){
		// find the unique contig
		Contig *uniqueContig = (*edgeIter)->GetOtherContig( repeat->GetOriginalRepeat() );
		int uniqueOri = (*edgeIter)->GetOrientationOfContig( uniqueContig );

		int repeatOri = (*edgeIter)->GetOrientationOfContig( repeat->GetOriginalRepeat() );
		if( repeatOri != repeat->GetOri() ){
			// reverse the orientation in the edge
			uniqueOri = GetOppositeOri( uniqueOri );
		}

		int repeatPos = (*edgeIter)->GetPositionOfContig( repeat );
		
		// check the edge
		int dis = 0;
		bool isSatisfied = false;
		int upperbound = (*edgeIter)->GetDis() + Configure::STD_TIMES * (*edgeIter)->GetStd();
		if( firstContig->m_contigName == repeat->GetName() && firstContig->m_ori == repeat->GetOri() 
		    && ( (repeatPos == START &&repeatOri == repeat->GetOri() ) || (repeatPos == END && repeatOri != repeat->GetOri() ) ) ){
			// check if the first contig can satisfy the edge
			for( list<SUBCONTIG*>::iterator subIter = region->begin(); subIter != region->end(); subIter++ ){
				// do not consider the first contig if it is the same as the current repeat
				if( subIter == region->begin() && (*subIter)->m_contigName == repeat->GetName() && (*subIter)->m_ori == repeat->GetOri() )
					continue;

				if( (*subIter)->m_contigName == uniqueContig->GetName() && (*subIter)->m_ori == uniqueOri ){
					// this edge can be satisfied
					isSatisfied = true;
					break;
				}
				
				dis += (*subIter)->m_length;
				if( dis >= upperbound ){
					// distance is longer than the upper-bound of this edge, so this edge cannot be satisfied
					break;
				}
			}
		}
		if( isSatisfied ){
			continue;
		}
		
		dis = 0;
		if( lastContig->m_contigName == repeat->GetName() && lastContig->m_ori == repeat->GetOri() 
		    && ( (repeatPos == END &&repeatOri == repeat->GetOri() ) || (repeatPos == START && repeatOri != repeat->GetOri() ) ) ){
			// check if the last contig can satisfy the edge
			for( list<SUBCONTIG*>::reverse_iterator subIter = region->rbegin(); subIter != region->rend(); subIter++ ){
				// do not consider the first contig if it is the same as the current repeat
				if( subIter == region->rbegin() && (*subIter)->m_contigName == repeat->GetName() && (*subIter)->m_ori == repeat->GetOri() )
					continue;

				if( (*subIter)->m_contigName == uniqueContig->GetName() && (*subIter)->m_ori == uniqueOri ){
					// this edge can be satisfied
					isSatisfied = true;
					break;
				}
				
				dis += (*subIter)->m_length;
				if( dis >= upperbound ){
					// distance is longer than the upper-bound of this edge, so this edge cannot be satisfied
					break;
				}
			}
		}
		if( !isSatisfied ){
			// found an edge which cannot be satisfied
			return false;
		}
	}

	// count the number of new concordant edges
	int numOfNewConcordantEdges = 0;
	// countNewConcordantEdges

	// all edges can be satisfied
	if( numOfNewDiscordantEdges == numOfNewConcordantEdges )
		return true;
	else
		return false;
}
*/
bool opera::CheckEdgeSatisfaction( list<PET*> *collectedEdges, list<SUBCONTIG*> *region, Contig *repeat ){
	// check after removing this repeat, if the number of discordant edges will increase

	// create a map of all subcontigs
	// map of forward order
	map<string, list<SUBCONTIG*>::iterator> mapOfContigs;
	// map of reverse order
	map<string, list<SUBCONTIG*>::reverse_iterator> reverseMapOfContigs;
	
	for( list<SUBCONTIG*>::iterator subIter = region->begin(); subIter != region->end(); subIter++ ){
		SUBCONTIG *currentContig = *subIter;
		// FIXME: save the status if the contig is repeat or not
		// fixed
		if( !m_graph->IsRepeat( currentContig->m_contigName ) ){
			// add unique contigs
			mapOfContigs.insert( pair<string, list<SUBCONTIG*>::iterator> (currentContig->m_contigName, subIter) );
		}
	}

	for( list<SUBCONTIG*>::reverse_iterator subIter = region->rbegin(); subIter != region->rend(); subIter++ ){
		SUBCONTIG *currentContig = *subIter;
		if( !m_graph->IsRepeat( currentContig->m_contigName ) ){
			// add unique contigs
			reverseMapOfContigs.insert( pair<string, list<SUBCONTIG*>::reverse_iterator> (currentContig->m_contigName, subIter) );
		}
	}

	// count the number of new discordant edges
	int numOfNewDiscordantEdges = 0;
	// for each edge, check if it is can be satisfied
	for( list<PET*>::iterator edgeIter = collectedEdges->begin(); edgeIter != collectedEdges->end(); edgeIter++ ){
		PET *currentEdge = *edgeIter;
		//cerr<<"checking edge: "<<currentEdge->GetOriString()<<endl;
		// FIXME: get the original contig (not the updated one)
		// FIXED
		Contig *uniqueContig = currentEdge->GetOtherContigOriginal( repeat->GetOriginalRepeat() );
		int uniqueOri = currentEdge->GetOrientationOfContigOriginal( uniqueContig );
		int uniquePos = currentEdge->GetPositionOfContigOriginal( uniqueContig );

		//int repeatOri = currentEdge->GetOrientationOfContigOriginal( repeat->GetOriginalRepeat() );

		SUBCONTIG *uniqueSubContig = *(mapOfContigs[ uniqueContig->GetName().c_str() ]); 

		// FIXME: get the original distance
		// fixed
		int upperbound = currentEdge->GetOriginalDistance() + Configure::STD_TIMES * currentEdge->GetStd();

		bool ifSatisfy = false;
		if( ( uniqueOri == uniqueSubContig->m_ori && uniquePos == START )
		    || ( uniqueOri != uniqueSubContig->m_ori && uniquePos == END ) ){
			// the repeat is on the right of the unique
			//cerr<<"right region\n";
			// find the unique contig in map
			list<SUBCONTIG*>::iterator contigIter = mapOfContigs[ uniqueSubContig->m_contigName ];
			contigIter++;
			int length = 0;
			
			while( contigIter != region->end() ){
				//cerr<<"contig: "<<(*contigIter)->m_contigName<<"\t"<<(*contigIter)->m_ori<<endl;
				if( (*contigIter)->m_contigName == repeat->GetName() && (*contigIter)->m_ori == repeat->GetOri() ){
					// this edge can be satisfied
					ifSatisfy = true;
					break;
				}

				length += (*contigIter)->m_length;
				if( length > upperbound + Configure::KMER )
					break;

				contigIter++;
			}
		}
		else{
			// the repeat is on the left of the unique
			//cerr<<"left region\n";
			// find the unique contig in reverse_map
			list<SUBCONTIG*>::reverse_iterator contigIter = reverseMapOfContigs[ uniqueSubContig->m_contigName ];
			contigIter++;
			int length = 0;
			
			while( contigIter != region->rend() ){
				if( (*contigIter)->m_contigName == repeat->GetName() && (*contigIter)->m_ori == repeat->GetOri() ){
					// this edge can be satisfied
					ifSatisfy = true;
					break;
				}

				length += (*contigIter)->m_length;
				if( length > upperbound + Configure::KMER )
					break;

				contigIter++;
			}
		}

		// if this edge can not be satisfied, add 1 to the numOfNewDiscordantEdges
		if( !ifSatisfy )
			numOfNewDiscordantEdges++;
	}

	// count the number of new concordant edges
	int numOfNewConcordantEdges = 0;
	// FIXME: create a list to save all the discordant edges
	list<PET*> allDiscordantEdges;
	for( list<PET*>::iterator edgeIter = allDiscordantEdges.begin(); edgeIter != allDiscordantEdges.end(); edgeIter++ ){
		// FIXME: get the original distance
		// FIXED
		int upperbound = (*edgeIter)->GetOriginalDistance() + Configure::STD_TIMES * (*edgeIter)->GetStd();

		// get the position of the unique contig
		Contig *uniqueContig, *otherContig;
		int uniqueContigOri, otherContigOri;
		int uniqueContigPos;
		//int otherContigPos;
		
		// FIXME: get the original contig
		// FIXED
		if( !(*edgeIter)->GetOriginalStartContig()->IsRepeat() ){
			uniqueContig = (*edgeIter)->GetOriginalStartContig();
			uniqueContigOri = (*edgeIter)->GetOriginalStartContigOri();
			uniqueContigPos = START;
			
			otherContig = (*edgeIter)->GetOriginalEndContig();
			otherContigOri = (*edgeIter)->GetOriginalEndContigOri();
			//otherContigPos = END;
		}
		else{
			uniqueContig = (*edgeIter)->GetOriginalEndContig();
			uniqueContigOri = (*edgeIter)->GetOriginalEndContigOri();
			uniqueContigPos = END;

			otherContig = (*edgeIter)->GetOriginalStartContig();
			otherContigOri = (*edgeIter)->GetOriginalStartContigOri();
			//otherContigPos = START;
		}

		SUBCONTIG *uniqueSubcontig = *(mapOfContigs[ uniqueContig->GetName() ]);
		bool ifSatisfy = false;

		if( ( uniqueContigOri == uniqueSubcontig->m_ori && uniqueContigPos == START )
		    || ( uniqueContigOri != uniqueSubcontig->m_ori && uniqueContigPos == END ) ){
			// if the other contig is on the right of uniqueContig, use mapOfContigs
			list<SUBCONTIG*>::iterator contigIter = mapOfContigs[ uniqueSubcontig->m_contigName ];
			contigIter++;
			int otherOri = otherContigOri;
			if( uniqueContigPos == END )
				otherOri = GetOppositeOri( otherContigOri );

			int length = 0;
			
			while( contigIter != region->end() ){
				if( (*contigIter)->m_contigName == otherContig->GetName() && (*contigIter)->m_ori == otherOri ){
					// this edge can be satisfied
					ifSatisfy = true;
					break;
				}

				length += (*contigIter)->m_length;
				if( length > upperbound + Configure::KMER )
					break;

				contigIter++;
			}
			
		}
		else{
			// else, use reverseMapOfContigs
			list<SUBCONTIG*>::reverse_iterator contigIter = reverseMapOfContigs[ uniqueSubcontig->m_contigName ];
			contigIter++;
			int otherOri = otherContigOri;
			if( uniqueContigPos == START )
				otherOri = GetOppositeOri( otherContigOri );

			int length = 0;
			
			while( contigIter != region->rend() ){
				if( (*contigIter)->m_contigName == otherContig->GetName() && (*contigIter)->m_ori == otherOri ){
					// this edge can be satisfied
					ifSatisfy = true;
					break;
				}

				length += (*contigIter)->m_length;
				if( length > upperbound + Configure::KMER )
					break;

				contigIter++;
			}
		}

		// if this edge can be satisfied, add 1 to the numOfNewConcordantEdges
		if( !ifSatisfy )
			numOfNewConcordantEdges++;
	}

	// all edges can be satisfied
	//cerr<<"new discordant edges: "<<numOfNewDiscordantEdges<<endl;
	//cerr<<"new concordant edges: "<<numOfNewConcordantEdges<<endl;

	if( numOfNewDiscordantEdges == numOfNewConcordantEdges ){
		// FIXME: label all the changes of the edges
		// labelEdgeChanges
		return true;
	}
	else
		return false;
}

// Given a repeat, collect the edges related with a set of contigs
// collectedEdges: all the selected edges
// possibleEdges: All the possible related edges
// repeat: the specific contig which is about to be removed
// region: the list of all possible related contigs
// direction: the direction of the edge from repeat (left or right)
void opera::CollectRelatedEdges( list<PET*> *collectedEdges, list<PET*> *possibleEdges, Contig *repeat, list<SUBCONTIG*> *region ){
	
	for( list<PET*>::iterator edgeIter = possibleEdges->begin(); edgeIter != possibleEdges->end(); edgeIter++ ){
		// do not consider the edge which is already used
		//if( (*edgeIter)->GetRepetitiveEdge() != NULL || (*edgeIter)->GetUniqueEdge() != NULL ){
		//	continue;
		//}

		// find the unique contig
		//cerr<<(*edgeIter)->GetOriString()<<endl;
		Contig *uniqueContig = (*edgeIter)->GetOtherContigOriginal( repeat->GetOriginalRepeat() );
		//cerr<<"unique contig is : "<<uniqueContig->GetName()<<endl;
		int uniqueOri = (*edgeIter)->GetOrientationOfContigOriginal( uniqueContig );
		//cerr<<"unique contig ori is: "<<uniqueOri<<endl;

		int repeatOri = (*edgeIter)->GetOrientationOfContig( repeat->GetOriginalRepeat() );
		if( repeatOri != repeat->GetOri() ){
			// reverse the orientation in the edge
			uniqueOri = GetOppositeOri( uniqueOri );
		}
		//cerr<<"repeat contig ori is: "<<uniqueOri<<endl;
		
		// find the contig in the region
		for( list<SUBCONTIG*>::iterator subIter = region->begin(); subIter != region->end(); subIter++ ){
			// check if the orientation is the same as in region
			if( (*subIter)->m_contigName == uniqueContig->GetName() && (*subIter)->m_ori == uniqueOri ){
				//cerr<<"subcontig is: "<<(*subIter)->m_contigName<<"\tori is: "<<(*subIter)->m_ori<<endl;
				// if it is the same, then add this edge to the colledtedEdges list
				collectedEdges->push_back( *edgeIter );
				break;
			}
		}
	}
}

void opera::CollectRelatedRightEdges( list<PET*> *collectedEdges, list<PET*> *possibleEdges, Contig *repeat, list<SUBCONTIG*> *region ){
	
	for( list<PET*>::iterator edgeIter = possibleEdges->begin(); edgeIter != possibleEdges->end(); edgeIter++ ){
		// do not consider the edge which is already used
		//if( (*edgeIter)->GetRepetitiveEdge() != NULL || (*edgeIter)->GetUniqueEdge() != NULL ){
		//	continue;
		//}
		if( (*edgeIter)->GetUniqueEdge() != NULL || (*edgeIter)->GetRepetitiveEdge() != NULL )
			continue;

		// find the unique contig
		//cerr<<(*edgeIter)->GetOriString()<<endl;
		Contig *uniqueContig = (*edgeIter)->GetOtherContigOriginal( repeat->GetOriginalRepeat() );
		//cerr<<"unique contig is : "<<uniqueContig->GetName()<<endl;
		int uniqueOri = (*edgeIter)->GetOrientationOfContigOriginal( uniqueContig );
		//cerr<<"unique contig ori is: "<<uniqueOri<<endl;

		int repeatOri = (*edgeIter)->GetOrientationOfContig( repeat->GetOriginalRepeat() );
		if( repeatOri != repeat->GetOri() ){
			// reverse the orientation in the edge
			uniqueOri = GetOppositeOri( uniqueOri );
		}
		//cerr<<"repeat contig ori is: "<<uniqueOri<<endl;
		
		// find the contig in the region
		int length = 0;
		for( list<SUBCONTIG*>::iterator subIter = region->begin(); subIter != region->end(); subIter++ ){
			// check if the orientation is the same as in region
			if( (*subIter)->m_contigName == uniqueContig->GetName() && (*subIter)->m_ori == uniqueOri ){
				//cerr<<"subcontig is: "<<(*subIter)->m_contigName<<"\tori is: "<<(*subIter)->m_ori<<endl;
				// if it is the same, then add this edge to the colledtedEdges list
				collectedEdges->push_back( *edgeIter );
				break;
			}

			length += (*subIter)->m_length;
			if( length > (*edgeIter)->GetOriginalDistance() + Configure::STD_TIMES * (*edgeIter)->GetStd() + Configure::KMER)
				break;
		}
	}
}

void opera::CollectRelatedLeftEdges( list<PET*> *collectedEdges, list<PET*> *possibleEdges, Contig *repeat, list<SUBCONTIG*> *region ){
	
	for( list<PET*>::iterator edgeIter = possibleEdges->begin(); edgeIter != possibleEdges->end(); edgeIter++ ){
		// do not consider the edge which is already used
		//if( (*edgeIter)->GetRepetitiveEdge() != NULL || (*edgeIter)->GetUniqueEdge() != NULL ){
		//	continue;
		//}

		// do not consider the edge which is already used
		if( (*edgeIter)->GetUniqueEdge() != NULL || (*edgeIter)->GetRepetitiveEdge() != NULL )
			continue;

		// find the unique contig
		//cerr<<(*edgeIter)->GetOriString()<<endl;
		Contig *uniqueContig = (*edgeIter)->GetOtherContigOriginal( repeat->GetOriginalRepeat() );
		//cerr<<"unique contig is : "<<uniqueContig->GetName()<<endl;
		int uniqueOri = (*edgeIter)->GetOrientationOfContigOriginal( uniqueContig );
		//cerr<<"unique contig ori is: "<<uniqueOri<<endl;

		int repeatOri = (*edgeIter)->GetOrientationOfContig( repeat->GetOriginalRepeat() );
		if( repeatOri != repeat->GetOri() ){
			// reverse the orientation in the edge
			uniqueOri = GetOppositeOri( uniqueOri );
		}
		//cerr<<"repeat contig ori is: "<<uniqueOri<<endl;
		
		// find the contig in the region
		int length = 0;
		for( list<SUBCONTIG*>::reverse_iterator subIter = region->rbegin(); subIter != region->rend(); subIter++ ){
			// check if the orientation is the same as in region
			if( (*subIter)->m_contigName == uniqueContig->GetName() && (*subIter)->m_ori == uniqueOri ){
				//cerr<<"subcontig is: "<<(*subIter)->m_contigName<<"\tori is: "<<(*subIter)->m_ori<<endl;
				// if it is the same, then add this edge to the colledtedEdges list
				collectedEdges->push_back( *edgeIter );
				break;
			}

			length += (*subIter)->m_length;
			if( length > (*edgeIter)->GetOriginalDistance() + Configure::STD_TIMES * (*edgeIter)->GetStd() + Configure::KMER )
				break;
		}
	}
}

// the pipeline of optimal scaffolding in opera
int main(int argc, char *argv[] )
{
	time_t start, end;
	struct timeval t_start,t_end;
	struct timeval t_startTemp,t_endTemp;
	time( &start );
#ifdef TIME
	gettimeofday( &t_start, NULL );
#endif
	
	// check parameters
	string configFileName = "";
	if( argc == 1 || argc == 3 || argc > 5 || string(argv[ 1 ]) == "-h" ){
		cout<<"Usage:\n";
		cout<<"  bin/opera <config-file>\n";
		cout<<"    OR\n";
		cout<<"  bin/opera <contig-file> <mapping-files> <output-folder>\n\n";
		cout<<"  <config-file> \tConfiguration file\n";
		cout<<"  <contig-file> \tMulti-fasta contig file\n";
		cout<<"  <mapping-files> \tComma-separated list of files containing mapping of paired-end reads\n";
		cout<<"  <output-folder> \tFolder to save the scaffold results\n";
		cout<<"  <samtools-dir> \tFolder which contains samtools binaries\n";
		cout<<"\n\nNOTE: Please refer to test_dataset/multiLib.config for detailed settings.\n";
		return -1;
	}
	else if( argc == 2 ){
		configFileName = string( argv[ 1 ] );
		
		// step 0: read config file
#ifdef TIME
		gettimeofday( &t_startTemp, NULL );
#endif
		cout<<"Step 1: Reading configuration file ..."<<endl;
		flush(cout);
		configureReader myConfig;
		if( myConfig.ReadConfigFile( configFileName ) == -1 ){
			cout<<"ERROR: The format of configuration file is not correct"<<endl; 
			return -1;
		}
#ifdef TIME
		gettimeofday( &t_endTemp, NULL );
		cout<<"Time Taken: "<<(t_endTemp.tv_sec - t_startTemp.tv_sec) + (t_endTemp.tv_usec - t_startTemp.tv_usec)/1000000.0<<" seconds"<<endl;
		flush(cout);
#endif
	}
	else if( argc == 5 || argc == 4 ){
		// get contig file name, mapping file name and output folder
		// step 0: read config file
#ifdef TIME
		gettimeofday( &t_startTemp, NULL );
#endif
		cout<<"Step 1: Setting parameters ..."<<endl;
		flush(cout);
		Configure::CONTIG_FILE = argv[ 1 ];

		// set mapping files
		vector<string> *tempName = new vector<string>;
		Split( argv[ 2 ], ",", tempName );
		for( int i = 0; i < (int) tempName->size(); i++ )
		{
			LibInfo *newLib = new LibInfo();
			newLib->SetFileName( tempName->at( i ) );
			
			if( IsSAM( newLib->GetFileName() ) )       
				newLib->SetMapType( SAM );
			Configure::MULTI_LIB_INFO->push_back( newLib );
		}
			
		Configure::OUTPUT_FOLDER = argv[ 3 ];
		mkdir (Configure::OUTPUT_FOLDER.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
		if( Configure::OUTPUT_FOLDER.substr( Configure::OUTPUT_FOLDER.length() - 1, 1 ) != Configure::DIRECTORY_SIGN )
			Configure::OUTPUT_FOLDER += Configure::DIRECTORY_SIGN;

		if( argc == 5){
                	Configure::SAMDIR = argv[ 4 ];
			Configure::SAMDIR = Configure::SAMDIR + "/";
		}else if( argc == 4){
			Configure::SAMDIR = "";
			//fprintf(stderr, "SAMDIR!!!!!!!!!!!! %s\n",Configure::SAMDIR.c_str());
		}		

#ifdef TIME
		gettimeofday( &t_endTemp, NULL );
		cout<<"Time Taken: "<<(t_endTemp.tv_sec - t_startTemp.tv_sec) + (t_endTemp.tv_usec - t_startTemp.tv_usec)/1000000.0<<" seconds"<<endl;
		flush(cout);
#endif
	}

	opera *m_opera = new opera();

	// check if all files exist
	if( !m_opera->CheckFileExist() ){
		return -1;
	}

	// step 0.5: check if the scaffold name will have conflicts
	m_opera->OpenLogFile();
	if( m_opera->CheckNameConfliction() ){
		// there is confliction
		cout<<"ERROR: There might be conflictions of scaffold name. Please specify another name using parameter: scaffold_name in configuraton file"<<endl;
		return -1;
	}

	// check contig file format
	//cerr<<"checking contig format\n";
	if( m_opera->CheckContigFormat() == -1 ){
		return -1;
	}

	// check the kmer size
	if( !Configure::USER_SPECIFY_KMER ){
		if( Configure::FILE_TYPE == VELVET ){
			if( !m_opera->GetKmer( Configure::CONTIG_FILE, VELVET ) ){
				// set the kmer to be 100
				//cout<<"Could not locate \"LastGraph\" file, the kmer size is set to be 100\n"; 
				Configure::KMER = 100;
			}
		}
		else if( Configure::FILE_TYPE == SOAP || Configure::FILE_TYPE == SOAP_CONTIG){
			if( !m_opera->GetKmer( Configure::CONTIG_FILE, SOAP ) ){
				// set the kmer to be 100
				//cout<<"Could not locate \".preGraphBasic\" file, the kmer size is set to be 100\n";
				Configure::KMER = 100;
			}
		}
		else{
			Configure::KMER = 100;
			//cout<<"The assembly format is neither Velvet nor SOAPdenovo. So the kmer size is set to be 100.\n";
		}
	}

	// output parameters to a file
	if( m_opera->PrintParameter( Configure::OUTPUT_FOLDER + "parameters" ) == -1 ){
		return -1;
	}

	// step 1: convert contig file
#ifdef TIME
	gettimeofday( &t_startTemp, NULL );
#endif
	cout<<"Step 2: Reading contig file ..."<<endl;
	flush(cout);
	ContigConverter myContigConverter;
	if( myContigConverter.ConvertContigFile( Configure::CONTIG_FILE, m_opera->m_graph, m_opera->m_libraries ) == -1 ){
		cout<<"ERROR: Converting contig file error!"<<endl;
		return -1;
	}
#ifdef TIME
	gettimeofday( &t_endTemp, NULL );
	cout<<"Time Taken: "<<(t_endTemp.tv_sec - t_startTemp.tv_sec) + (t_endTemp.tv_usec - t_startTemp.tv_usec)/1000000.0<<" seconds"<<endl;
	flush(cout);
#endif

	if( Configure::REAL_POSITION_FILE != "" ){
		// read the real position file of all the unique contigs
		if( m_opera->ReadRealPositionFile( Configure::REAL_POSITION_FILE ) == -1 ){
			return -1;
		}
	}

	// step 2: convert mapping file
#ifdef TIME
	gettimeofday( &t_startTemp, NULL );
#endif
	cout<<"Step 3: Reading mapping file ..."<<endl;
	flush(cout);
	MapConverter myMapConverter( m_opera->m_graph );
	/*if( myMapConverter.Analyze( Configure::MAP_FILE ) == -1 ){
		cout<<"ERROR: Converting mapping file error!"<<endl;
		return -1;
	}*/
	if( myMapConverter.AnalyzeMultiLib( Configure::MAP_FILE, m_opera->m_libraries ) == -1 ){
		cout<<"ERROR: Converting mapping file error!"<<endl;
		return -1;
	}

#ifdef TIME
	gettimeofday( &t_endTemp, NULL );
	cout<<"Time Taken: "<<(t_endTemp.tv_sec - t_startTemp.tv_sec) + (t_endTemp.tv_usec - t_startTemp.tv_usec)/1000000.0<<" seconds"<<endl;
	flush(cout);
#endif

	////////////////////////////////////
	// Test table construction of gap sizes
	/*
	for( int contigLengthSum = 1000; contigLengthSum <= 16000; contigLengthSum += 1000 ){
		cout<<"Constructing gap table for contig length: "<<contigLengthSum<<"...\n";
		GapCorrecter *gc = new GapCorrecter( 16000, contigLengthSum, 10000, 1000 );
		gc->CalculateTable();
		cout<<"Printing gap table...\n";
		gc->PrintTable( "gapTable" + itos( contigLengthSum) );
	}
	return -1;
	*/
	///////////////////////////////////

	// step 3: bundle
#ifdef TIME
	gettimeofday( &t_startTemp, NULL );
#endif
	cout<<"Step 4: Bundling paired-end reads ..."<<endl;
	flush(cout);
	/*if( Configure::MAP_TYPE != OPERA ){
		if( myMapConverter.Bundle() == -1 ){
			cout<<"ERROR: Bundling error!"<<endl;
			return -1;
		}
	}*/

	if( myMapConverter.BundleMultiLib( m_opera->m_libraries, m_opera->m_discardEdgeFile, m_opera->m_superLibrary,
					   m_opera->m_conflictingEdges ) == -1 ){
		cout<<"ERROR: Bundling error!"<<endl;
		return -1;
	}

	// output conflicting edges
	if( m_opera->OutputConflictingEdges( Configure::OUTPUT_FOLDER + "conflictingEdges" ) == -1 )
		return -1;

#ifdef TIME
	gettimeofday( &t_endTemp, NULL );
	cout<<"Time Taken: "<<(t_endTemp.tv_sec - t_startTemp.tv_sec) + (t_endTemp.tv_usec - t_startTemp.tv_usec)/1000000.0<<" seconds"<<endl;
	flush(cout);
#endif

	list<PetLibrary*>::iterator libIter;
	list<PetLibrary*> *libraries;

	if( !Configure::SIMULTANEOUSLY_HANDLE_LIBRARIES ){
		libIter = m_opera->m_libraries->begin();
		libraries = m_opera->m_libraries;
	}
	else{
		libIter = m_opera->m_superLibrary->begin();
		libraries = m_opera->m_superLibrary;
	}
	if( !Configure::ONLY_BUNDLE )
	{
		// traverse each library
		int libNum = 0;
	
		// step 4: scaffolder
#ifdef TIME
		gettimeofday( &t_startTemp, NULL );
#endif
		cout<<"Step 5: Scaffolding ..."<<endl;
		flush(cout);

		while( libIter != libraries->end() ){
			m_opera->m_totalNumberOfSubgraph = 0;
			m_opera->m_totalNumberOfIncreasedSubgraph = 0;
			m_opera->m_totalNumberOfNodes = 0;
			m_opera->m_totalIncreasedNodes.clear();

			libNum++;
			cout<<"Scaffolding with library: "<<(*libIter)->GetFileName()<<endl;

			PetLibrary *currentLib = *libIter;
			Configure::LIB_MEAN = currentLib->GetMean();
			Configure::LIB_STD = currentLib->GetStd();
			Configure::UPPERBOUND = Configure::LIB_MEAN + Configure::STD_TIMES * Configure::LIB_STD;

			// initialize opera for a new run
			m_opera->Initialize();

			// move the edges to graph
			m_opera->MoveEdgesToGraph( currentLib );
			//cout<<"finish moving edges to graph\n";

			// delete current library
			delete currentLib;

			// split the contigs according to their types
			m_opera->m_graph->InitializeContig( libNum );
			m_opera->m_graph->GenerateBorderAndScaffold();
			
#ifdef DEBUG
			cout<<"finish generating border contig\n";
#endif

			//m_opera->m_graph->CountSubgraph();
			//return -1;

			if( m_opera->StartScaffold() == -1 ){
				cout<<"ERROR: Scaffolding error!"<<endl;
				return -1;
			}

			// check the percentage of subgraphs increasing threshold
			double percentageOfIncreasedSubgraph = (double) m_opera->m_totalNumberOfIncreasedSubgraph / (double) m_opera->m_totalNumberOfSubgraph;
			cout<<"Percentage of subgraphs which need to increase the cluster size threhold: "<<percentageOfIncreasedSubgraph<<endl;
			if( percentageOfIncreasedSubgraph > Configure::PERCENTAGE_OF_INCREASING_THRESHOLD_GRAPHS ){
				cout<<"WARNING: more than "<<Configure::PERCENTAGE_OF_INCREASING_THRESHOLD_GRAPHS<<"% of subgraphs need to increase cluster size threshold locally. It is recommended to specify a bigger global threshold in configuration file using \"cluster_size=\""<<endl;  
			}

			double percentageOfIncreasedNode = (double) m_opera->m_totalIncreasedNodes.size() / (double) m_opera->m_totalNumberOfNodes;
			cout<<"Percentage of nodes which need to increase the cluster size threhold: "<<percentageOfIncreasedNode<<endl;
			if( percentageOfIncreasedNode > Configure::PERCENTAGE_OF_INCREASING_THRESHOLD_GRAPHS ){
				cout<<"WARNING: more than "<<Configure::PERCENTAGE_OF_INCREASING_THRESHOLD_GRAPHS<<"% of nodes need to increase cluster size threshold locally. It is recommended to specify a bigger global threshold in configuration file using \"cluster_size=\""<<endl;  
			}


#ifdef TIME
			gettimeofday( &t_endTemp, NULL );
			cout<<"Time Taken: "<<(t_endTemp.tv_sec - t_startTemp.tv_sec) + (t_endTemp.tv_usec - t_startTemp.tv_usec)/1000000.0<<" seconds"<<endl;
			flush(cout);
#endif

			// remove first library
			libIter = libraries->erase( libIter );
		}

	// step 6: finalize results
#ifdef TIME
		gettimeofday( &t_startTemp, NULL );
#endif
		cout<<"Step 6: Outputing results ..."<<endl;
		flush(cout);
		if( m_opera->OutputScaffold( Configure::OUTPUT_FOLDER + "scaffolds.scaf" ) == -1 ){
			cout<<"ERROR: Output scaffolds result error!"<<endl;
			return -1;
		}
#ifdef TIME
		//cout<<"finish output .scaf file. Now generating sequence file..."<<endl;
		//flush(cout);
#endif

		//cerr<<"done outputting scaffolds\n";
		m_opera->SortScaffold();

		if( m_opera->OutputUnhappyEdges( Configure::OUTPUT_FOLDER + "discordantEdges" ) == -1 ){
			cout<<"ERROR: Output discordant edges error!"<<endl;
			return -1;
		}
#ifdef TIME
		gettimeofday( &t_endTemp, NULL );
		cout<<"Time Taken: "<<(t_endTemp.tv_sec - t_startTemp.tv_sec) + (t_endTemp.tv_usec - t_startTemp.tv_usec)/1000000.0<<" seconds"<<endl;
		gettimeofday( &t_end, NULL);
#endif

		// print the statistics
		cout<<m_opera->GetStats()<<endl;
	} // end of if( !only bundle )

	bool printTime = false;
	time( &end );
	cout<<"Total running time is ";
	int days = 0, hours = 0, minutes = 0, seconds = 0;
	double dif = difftime( end, start );
	if( dif > 60*60*24 ){
		days = (int)( dif / (60*60*24) );
		dif -= days * (60*60*24);
	}
	if( dif > 60*60 ){
		hours = (int) (dif / (60*60));
		dif -= hours * (60 * 60);
	}
	if( dif > 60 ){
		minutes = (int) (dif / 60);
		dif -= minutes * 60;
	}

	seconds = (int)dif;
	if( days != 0 ){
		printTime = true;
		cout<<days<<" days ";
	}
	if( hours != 0 ){
		printTime = true;
		cout<<hours<<" hours ";
	}
	if( minutes != 0 ){
		printTime = true;
		cout<<minutes<<" minutes ";
	}
	if( seconds != 0 ){
		printTime = true;
		cout<<seconds<<" seconds ";
	}
	if( !printTime ){
		cout<<"0 seconds";
	}
	cout<<endl;
	
	//cout<<"Percentage of top unassigned nodes: "<<((double)m_opera->m_numOfTopUnassignedNodes/(double)m_opera->m_totalPartialScaffold);
	//cout<<"\t("<<m_opera->m_numOfTopUnassignedNodes<<"/"<<m_opera->m_totalPartialScaffold<<")\n";

	delete m_opera;
	
	cout<<"The results are in "<<Configure::OUTPUT_FOLDER.substr( 0, Configure::OUTPUT_FOLDER.length() - 1 )<<endl;

	Configure::destroy();

#ifdef DEBUG
	system("pause");
#endif
}

