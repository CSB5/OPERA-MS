#include "Graph.h"

Graph::Graph(void)
{
	//m_contigNameMap = new map<string, int>;
	m_contigNameHashMap = new hash_map<const char*, int, hash<const char*>, eqName>;
	m_numberOfContigs = 0;
	m_numberOfEdges = 0;

	m_contigsList = new list<Contig*>;		
	m_scaffoldsList = new list<Contig*>;
	m_borderContigsList = new list<Contig*>;

	ifInitializeContig = false;

	m_contigID = new map<int, Contig*>();
	m_contigNames = new set<string>;

	m_repeatContigSet = new hash_map<const char*, bool, hash<const char*>, eqName>;
}

Graph::~Graph(void)
{
	if( ifInitializeContig )
		delete[] m_contigsArray;
	DeleteContigs( m_contigsList );
	delete m_contigsList;
	DeleteContigs( m_scaffoldsList );
	delete m_scaffoldsList;
	m_borderContigsList->clear();
	delete m_borderContigsList;

	//m_contigNameMap->clear();
	//delete m_contigNameMap;
	m_contigNameHashMap->clear();
	delete m_contigNameHashMap;

	m_contigNames->clear();
	delete m_contigNames;

	m_contigID->clear();
	delete m_contigID;

	m_repeatContigSet->clear();
	delete m_repeatContigSet;
}

// save contig name
void Graph::SaveContigName( string name ){
	m_contigNames->insert( name );
}

// check if a contig name exist
bool Graph::IfExistInContig( string name ){
	if( m_contigNames->find( name ) == m_contigNames->end() )
		return false;
	else
		return true;
}

void Graph::DeleteContigs( list<Contig*> *contigs ){
	list<Contig*>::iterator iter = contigs->begin();
	while( iter != contigs->end() ){
		Contig* temp = *iter;
		iter = contigs->erase( iter );
		delete temp;
	}
}

// set total number of contigs
void Graph::SetContigNum( int n ){
	m_contigsArray = new Contig*[ n ];
	ifInitializeContig = true;
}

// Add non-repeat contigs
void Graph::AddContig( Contig *c ){
	// add to array
	m_contigsArray[ m_numberOfContigs ] = c;
	c->SetID( m_numberOfContigs );
	//m_contigNameMap->insert( pair<string, int>(c->GetName(), c->GetID() ) );
	m_contigNameHashMap->insert( pair<const char*, int>(c->GetName().c_str(), c->GetID() ) );
	m_numberOfContigs++;

	// add to list
	m_contigsList->push_front( c );
	c->SetListPos( m_contigsList->begin() );
}

// get the contig in position pos
Contig* Graph::GetContigUsingPos( int pos )
{
	return m_contigsArray[ pos ];
}

// initialize contig list for a new library
void Graph::InitializeContig( int libNum ){
	// do not need to initialize for the first library
	if( libNum == 1 )
		return;

	m_contigNameHashMap->clear();

	m_numberOfContigs = 0;
	list<Contig*>::iterator scaIter = this->m_scaffoldsList->begin();
	//cout<<"number of contigs: "<<m_scaffoldsList->size()<<endl;
	while( scaIter != m_scaffoldsList->end() ){
		this->m_contigsList->push_front( *scaIter );
		(*scaIter)->SetListPos( m_contigsList->begin() );
		(*scaIter)->Initialize();
		m_contigsArray[ m_numberOfContigs ] = *scaIter;
		(*scaIter)->SetID( m_numberOfContigs );
		m_contigNameHashMap->insert( pair<const char*, int>((*scaIter)->GetName().c_str(), (*scaIter)->GetID() ) );
		
		m_numberOfContigs++;
		scaIter = m_scaffoldsList->erase( scaIter );
	}
}

// output the contigs into a new file
// return -1 if failed
int Graph::OutputContigs( string fileName ){
	ofstream contigWriter( fileName.c_str() );

	if( contigWriter == NULL )
	{
		cout<<"ERROR: Cannot open "<<fileName<<" file"<<endl;
		return -1;
	}

	string head = "Contig ID\tLength\tCoverage\n";
	contigWriter.write( head.c_str(), head.length() );
	
	list<Contig*>::const_iterator iter;
	for( iter = m_contigsList->begin(); iter != m_contigsList->end(); iter++ ){
		string line = (*iter)->GetName() + "\t";
		char buffer[30];
		sprintf( buffer, "%.0f\t%.0f\n", (*iter)->GetLength(), (*iter)->GetCov() );
		string tempString( buffer );
		line += tempString;
		//line += ( itos( (*iter)->GetLength() )+ "\n");
		contigWriter.write( line.c_str(), line.length() );
	}

	return 1;
}

// get contig index according to contig name
// return -1 if not find
int Graph::GetContigIndex( string contigName ){
	//map<string, int>::const_iterator pos = m_contigNameMap->find( contigName );
	hash_map<const char*, int, hash<const char*>, eqName>::iterator pos = m_contigNameHashMap->find( contigName.c_str() );
	if( pos == m_contigNameHashMap->end() )
		return -1;
	else
		return pos->second;

	
}

// get contig in specific index
Contig * Graph::GetContig( int pos ){
	return m_contigsArray[ pos ];
}

// add edge to graph
void Graph::AddEdge( PET *p ){
	p->SetID( m_numberOfEdges++ );
	p->GetStartContig()->AddEdge( p );
	p->GetEndContig()->AddEdge( p );
}

// add edges from all libraries to graph
void Graph::AddEdgeMultiLib( PET *p ){
	//p->SetID( m_numberOfEdges++ );
	p->GetStartContig()->AddEdgeMultiLib( p );
	p->GetEndContig()->AddEdgeMultiLib( p );
}

// generate border contig and singletons
void Graph::GenerateBorderAndScaffold(){
	//cerr<<"generating border contig\n";
	int threshold = Configure::UPPERBOUND;
	//cerr<<"uppderbound is "<<threshold<<"\n";

	list<Contig*>::iterator iter = m_contigsList->begin();
	//cerr<<"First contig is :"<<(*iter)->GetName()<<endl;
	while( iter != m_contigsList->end() ){
		//cerr<<(*iter)->GetName()<<"\t"<<(*iter)->GetLength()<<"\ttype: "<<(*iter)->GetContigType()<<endl;
		if( (*iter)->IsSingleton() || (*iter)->IsRepeat() ){
			// save to scaffold list (all repeats are singletons)
			m_scaffoldsList->push_back( *iter );
			iter = m_contigsList->erase( iter );
			//cerr<<"is singleton\n";
		}
		else if( (*iter)->GetLength() >= threshold && !(*iter)->IsRepeat() ){
			// save to border contig List
			//cerr<<"border contig: "<<(*iter)->GetName()<<endl;
			m_borderContigsList->push_front( *iter );
			(*iter)->SetBorderListPos( m_borderContigsList->begin() );
			iter++;
			//cerr<<"is border contig\n";
		}
		else
			iter++;
	}

	//cout<<"there are "<<this->m_contigsList->size()<<" contigs with edge\n";
}

// update the graph using results of subgraph
// handle repeats: repeats should not be output as singletons; edge to repeat should be updated rather than deleted.
//                 repeats should not be used as "super contigs"
void Graph::UpdateGraph( list<Contig*> *subgraph, ScaffoldResult **scaffold, FILE* discardEdgeFile ){
	bool isSpecial = false;
	
#ifdef DEBUG
	//cerr<<"start update graph\n";
	//	cerr<<"handle each contig\n";
#endif
	// deal with each contigs
	list<Contig*> *tempListOfContig = new list<Contig*>;	// save all contigs needed to be deleted in subgraph

	list<Contig*>::iterator contigIter = subgraph->begin();
	while( contigIter != subgraph->end() ){
		//cout<<"scaffold\n"<<scaffold[ (*contigIter)->GetScaffoldID() ]->GetScaffoldString()<<endl;
		//	cerr<<"contig:\n"<<(*contigIter)->GetScaffoldString()<<endl;
		if( isSpecial )
			cerr<<"Handle contig: "<<(*contigIter)->GetName()<<endl;

		bool isFirstContig = false;
		bool sameOri = false;
		// check if need to reverse left and right edges
		sameOri = SameOri( *contigIter, scaffold[ (*contigIter)->GetScaffoldID() ]->GetScaffoldString() );

		//if( scaffold[ (*contigIter)->GetScaffoldID() ]->GetID() == -1 ){
		if( scaffold[ (*contigIter)->GetScaffoldID() ]->GetID() == -1 
		    || m_contigsArray[ scaffold[ (*contigIter)->GetScaffoldID() ]->GetID() ]->IsRepeat() ){
			// do not use a repeat as the combined nodes, unless it is the only node

			if( !sameOri && !(*contigIter)->IsRepeat()  )
				(*contigIter)->ReverseEdges();
			//if( IfSameOriInScaffoldAndEdge( *contigIter, MINUS, scaffold[ (*contigIter)->GetScaffoldID() ] ) )
			//	(*contigIter)->ReverseEdges();
			isFirstContig = true;
			scaffold[ (*contigIter)->GetScaffoldID() ]->SetID( (*contigIter)->GetID() );
			if( !(*contigIter)->IsRepeat() ){
				(*contigIter)->SetScaffoldString( scaffold[ (*contigIter)->GetScaffoldID() ]->GetScaffoldString() );
				(*contigIter)->SetLength( scaffold[ (*contigIter)->GetScaffoldID() ]->GetLength() );
				(*contigIter)->SetCov( scaffold[ (*contigIter)->GetScaffoldID() ]->GetCov() );
				(*contigIter)->SetLengthWithGap( scaffold[ (*contigIter)->GetScaffoldID() ]->GetLengthWithGap() );
			}
		}

		//cout<<(*contigIter)->GetScaffoldString()<<endl;

		//if( isSpecial )
		//	cerr<<"Deleting edge of current library\n";

		if( !(*contigIter)->IsRepeat() ){
#ifdef DEBUG
			//cerr<<"delete edges\n";
#endif

			// traverse each edge
			DeleteEdges( *contigIter, (*contigIter)->GetLeftEdges(),
				     isFirstContig, scaffold[ (*contigIter)->GetScaffoldID() ], sameOri );

			DeleteEdges( *contigIter, (*contigIter)->GetRightEdges(),
				     isFirstContig, scaffold[ (*contigIter)->GetScaffoldID() ], sameOri );
		

#ifdef DEBUG
			//cerr<<"Deleting edges of other libraries\n";
#endif

			// for multiple libraries
			// update the edges in other libraries
			DeleteEdgesMultiLib( *contigIter, (*contigIter)->GetLeftEdgeMultiLib(),
					     isFirstContig, scaffold[ (*contigIter)->GetScaffoldID() ], sameOri, discardEdgeFile );
			
			DeleteEdgesMultiLib( *contigIter, (*contigIter)->GetRightEdgeMultiLib(),
					     isFirstContig, scaffold[ (*contigIter)->GetScaffoldID() ], sameOri, discardEdgeFile );
		}

#ifdef DEBUG
		//cerr<<"delete current contig\n";
#endif

		// delete contig
		if( !(*contigIter)->IsRepeat() ){
			if( isFirstContig ){
				// update border contig pointer
				if( !(*contigIter)->IsBorderContig() ){
					m_borderContigsList->push_front( *contigIter );
					(*contigIter)->SetBorderListPos( m_borderContigsList->begin() );
				}
				contigIter++;
			}
			else{
				// delete current contig and remove the pointer
				m_contigsList->erase( (*contigIter)->GetListPos() );
				if( (*contigIter)->IsBorderContig() )
					m_borderContigsList->erase( (*contigIter)->GetBorderListPos() );
				m_contigsArray[ (*contigIter)->GetID() ] = NULL;
				//Contig* temp = *contigIter;
				tempListOfContig->push_back( *contigIter );
				contigIter = subgraph->erase( contigIter );
				//delete temp;
			}
		}
		else{
			contigIter++;
		}
	}


#ifdef DEBUG
	//cerr<<"deleting contigs\n";
#endif

	// delete contigs
	contigIter = tempListOfContig->begin();
	while( contigIter != tempListOfContig->end() ){
		Contig *temp = *contigIter;
		delete temp;
		contigIter = tempListOfContig->erase( contigIter );
	}
	delete tempListOfContig;
	if( isSpecial )
		cerr<<"finish deleting contigs\n";

	
#ifdef DEBUG
	//cerr<<" check remaining new contigs if they are singletons\n";
#endif
	// check remaining new contigs if they are singletons
	contigIter = subgraph->begin();
	while( contigIter != subgraph->end() ){
		(*contigIter)->SetNotInSubgraph();
		/*if( !(*contigIter)->IsRepeat() ){
			// FIXME: check edges to see if two repeats are adjacent to each other
			// fixed
			CheckEdgesForAdjacentRepeat( (*contigIter)->GetLeftEdges() );
			CheckEdgesForAdjacentRepeat( (*contigIter)->GetRightEdges() );
			CheckEdgesForAdjacentRepeat( (*contigIter)->GetLeftEdgeMultiLib() );
			CheckEdgesForAdjacentRepeat( (*contigIter)->GetRightEdgeMultiLib() );
		}
		*/

		if( !(*contigIter)->HasEdge() && !(*contigIter)->IsRepeat() ){
			// all the non repetitive contigs without any edges are removed to singleton list
			m_scaffoldsList->push_back( *contigIter );
			m_contigsArray[ (*contigIter)->GetID() ] = NULL;
			if( (*contigIter)->IsBorderContig() )
				m_borderContigsList->erase( (*contigIter)->GetBorderListPos() );
			m_contigsList->erase( (*contigIter)->GetListPos() );
		}
		/*else if( (*contigIter)->HasEdge() && !(*contigIter)->IsRepeat() ){
			// FIXME: check edges to see if two repeats are adjacent to each other
			// fixed
			CheckEdgesForAdjacentRepeat( (*contigIter)->GetLeftEdges() );
			CheckEdgesForAdjacentRepeat( (*contigIter)->GetRightEdges() );
			CheckEdgesForAdjacentRepeat( (*contigIter)->GetLeftEdgeMultiLib() );
			CheckEdgesForAdjacentRepeat( (*contigIter)->GetRightEdgeMultiLib() );
		}
		*/

		contigIter++;
	}
	if( isSpecial )
		cerr<<"finish checking singleton contigs\n";
}

// check if two contigs contain the same repeats near each other, if yes, modify the distance
void Graph::CheckEdgesForAdjacentRepeat( list<PET*> *edges ){
	for( list<PET*>::iterator edgeIter = edges->begin(); edgeIter != edges->end(); edgeIter++ ){
		PET *currentEdge = *edgeIter;
		if( currentEdge->GetStartContig()->IsRepeat() || currentEdge->GetEndContig()->IsRepeat() ){
			// only check the edges between unique contigs
			continue;
		}

		SUBCONTIG *firstContig, *secondContig;
		if( currentEdge->GetStartContigOri() == PLUS ){
			firstContig = currentEdge->GetStartContig()->GetEndSubContig( END, false );
		}
		else{
			firstContig = currentEdge->GetStartContig()->GetEndSubContig( START, true );
		}

		if( currentEdge->GetEndContigOri() == PLUS ){
			secondContig = currentEdge->GetEndContig()->GetEndSubContig( START, false );
		}
		else{
			secondContig = currentEdge->GetEndContig()->GetEndSubContig( END, true );
		}

		int oriDis = currentEdge->GetDis();
		if( firstContig->m_contigName == secondContig->m_contigName 
		    && firstContig->m_ori == secondContig->m_ori ){
			// increase the distance of this edge
			currentEdge->SetDis( currentEdge->GetDis() + firstContig->m_length );
		}

		if( currentEdge->GetDis() > 30000 ){
			cerr<<"Error: \n";
			cerr<<currentEdge->GetOriString()<<endl;
			cerr<<"Ori distance: "<<oriDis<<"\tnew distance: "<<currentEdge->GetDis()<<endl;
			cerr<<"firstContig: "<<firstContig->m_contigName<<"\tlength: "<<firstContig->m_length<<endl;
			cerr<<"start contig: "<<currentEdge->GetStartContig()->GetName()<<" "<<currentEdge->GetStartContigOri()<<endl;
			cerr<<currentEdge->GetStartContig()->GetScaffoldString()<<endl;
			cerr<<"end contig: "<<currentEdge->GetEndContig()->GetName()<<" "<<currentEdge->GetEndContigOri()<<endl;
			cerr<<currentEdge->GetEndContig()->GetScaffoldString()<<endl;
		}

		delete firstContig;
		delete secondContig;
		
	}
}

void Graph::CheckEdgesForAdjacentRepeat( PET *currentEdge ){
	if( currentEdge->GetStartContig()->IsRepeat() || currentEdge->GetEndContig()->IsRepeat() ){
		// only check the edges between unique contigs
		return;
	}

	SUBCONTIG *firstContig, *secondContig;
	if( currentEdge->GetStartContigOri() == PLUS ){
		firstContig = currentEdge->GetStartContig()->GetEndSubContig( END, false );
	}
	else{
		firstContig = currentEdge->GetStartContig()->GetEndSubContig( START, true );
	}
	
	if( currentEdge->GetEndContigOri() == PLUS ){
		secondContig = currentEdge->GetEndContig()->GetEndSubContig( START, false );
	}
	else{
		secondContig = currentEdge->GetEndContig()->GetEndSubContig( END, true );
	}
	
	int oriDis = currentEdge->GetDis();
	if( firstContig->m_contigName == secondContig->m_contigName 
	    && firstContig->m_ori == secondContig->m_ori ){
		// increase the distance of this edge
		currentEdge->SetDis( currentEdge->GetDis() + firstContig->m_length );
	}
	
	if( currentEdge->GetDis() > 30000 ){
		cerr<<"Error: \n";
		cerr<<currentEdge->GetOriString()<<endl;
		cerr<<"Ori distance: "<<oriDis<<"\tnew distance: "<<currentEdge->GetDis()<<endl;
		cerr<<"firstContig: "<<firstContig->m_contigName<<"\tlength: "<<firstContig->m_length<<endl;
		cerr<<"start contig: "<<currentEdge->GetStartContig()->GetName()<<" "<<currentEdge->GetStartContigOri()<<endl;
		cerr<<currentEdge->GetStartContig()->GetScaffoldString()<<endl;
		cerr<<"end contig: "<<currentEdge->GetEndContig()->GetName()<<" "<<currentEdge->GetEndContigOri()<<endl;
		cerr<<currentEdge->GetEndContig()->GetScaffoldString()<<endl;
	}

	delete firstContig;
	delete secondContig;
}

// check if contig c has the same orientation as in scaffold
bool Graph::SameOri( Contig *c, string scaffold ){
	// get old contig information
	vector<string> lineOfScaf;
	Split( c->GetScaffoldString(), "\n", &lineOfScaf );
	list<string> line;
	Split( c->GetScaffoldString(), "\t", &line );

	//cerr<<c->GetScaffoldString()<<endl;
	//cerr<<lineOfScaf.size()<<endl;
	list<string>::iterator iter = line.begin();
	
	if( !c->IsRepeat() ){
		int i = 0;
		while( IsRepeat( *iter ) ){
			i++;
			//cerr<<"current line: "<<lineOfScaf.at( i )<<endl;
			//line.clear();
			Split( lineOfScaf.at( i ), "\t", &line );
			//cerr<<line.front()<<endl;
			iter = line.begin();
			//cerr<<"i: "<<i<<endl;
			//cerr<<"contig: "<<*iter<<endl;
		}
	}
	
	string oldContigName = *iter;
	iter++;
	string oldContigOri = *iter;
		
	/*if( oldContigName == "1593" )
	{
		cout<<"contig:\n"<<c->GetScaffoldString()<<endl;
		cout<<"scaffold:\n"<<scaffold<<endl;
	}*/

	// get new contig information
	list<string> lines;
	Split( scaffold, "\n", &lines );

	iter = lines.begin();
	while( iter != lines.end() ){
		vector<string> content;
		Split( *iter, "\t", &content );
		if( content.at( 0 ) == oldContigName ){
			if( content.at( 1 ) == oldContigOri )
				return true;
			else
				return false;
		}
		iter++;
	}

	return false;
}

// release the memory of edges
void Graph::DeleteEdges( Contig *c, list<PET*> *edges, bool isFirst, ScaffoldResult *scaffold, bool ifSame ){
	list<PET*>::iterator iter = edges->begin();
	while( iter != edges->end() ){
		// if this edge is within the subgraph
		if( (*iter)->IfInSubgraph() ){
			if( !(*iter)->IfVisited() ){
				Contig *otherContig = (*iter)->GetOtherContig( c );
				int otherPos = (*iter)->GetPositionOfContig( otherContig );
				int otherOri = (*iter)->GetOrientationOfContig( otherContig );
				if( !otherContig->IsRepeat() ){
					// if it is the first time, update the status of edge
					(*iter)->VisitEdge();
					iter = edges->erase( iter );
				}
				else{
					// directly delete this edge
					PET *temp = *iter;
					iter = edges->erase( iter );

					// delete from original repeat
					if( (otherOri == PLUS && otherPos == START) || (otherOri == MINUS && otherPos == END) ){
						// right edge
						if( temp->GetUniqueEdge() != NULL ){
							// if it is a unique edge
							temp->GetUniqueEdge()->SetOriginalEdge( NULL );
							// delete unique edge?
							// delete temp->GetUniqueEdge();

							temp->SetUniqueEdge( NULL );	
						}
						else if( temp->GetRepetitiveEdge() != NULL ){
							// if it is has a repeat edge

							//temp->GetRepetitiveEdge()->SetOriginalEdge( NULL );
							// delete repetitive edge
							PET *repetitiveEdge = temp->GetRepetitiveEdge();
							Contig *copyOfRepetitiveContig = repetitiveEdge->GetOtherContig( c );
							copyOfRepetitiveContig->RemoveOriginalRepetitiveEdge( temp, RIGHT );
							temp->GetRepetitiveEdge()->SetOriginalEdge( NULL );
							temp->SetRepetitiveEdge( NULL );
						}

						otherContig->RemoveOriginalRepetitiveEdge( temp, RIGHT );
					}
					else{
						// left edge
						if( temp->GetUniqueEdge() != NULL ){
							// if it is a unique edge
							temp->GetUniqueEdge()->SetOriginalEdge( NULL );
							// delete unique edge?
							// delete temp->GetUniqueEdge();

							temp->SetUniqueEdge( NULL );
						}
						else if( temp->GetRepetitiveEdge() != NULL ){
							// if it is has a repeat edge
							//temp->GetRepetitiveEdge()->SetOriginalEdge( NULL );
							// delete repetitive edge
							PET *repetitiveEdge = temp->GetRepetitiveEdge();
							Contig *copyOfRepetitiveContig = repetitiveEdge->GetOtherContig( c );
							copyOfRepetitiveContig->RemoveOriginalRepetitiveEdge( temp, LEFT );
							temp->GetRepetitiveEdge()->SetOriginalEdge( NULL );
							temp->SetRepetitiveEdge( NULL );
						}

						otherContig->RemoveOriginalRepetitiveEdge( temp, LEFT );
					}

					// fix me: memory leak
					delete temp;
				}
			}
			else{
				// the second time, delete
				PET *temp = *iter;
				iter = edges->erase( iter );
				delete temp;
			}
		}
		else if( c->GetScaffoldID() != -1){
			Contig *tempOtherContig = (*iter)->GetOtherContig( c );
			//cout<<"Edge is needed to be updated"<<endl;
			//cout<<"before: "<<(*iter)->ToString()<<endl;
			//cout<<"if same: "<<ifSame<<endl;
			//cout<<c->GetName()<<endl;
			// if this edge connects outside the subgraph, update it
			if( isFirst ){
				// meet the first contig in a scaffold, replace it
				if( scaffold->GetID() == -1 ){
					scaffold->SetID( c->GetID() );
					// update length of contig
					c->SetLength( scaffold->GetLength() );
				}
			}

			// replace edges
			(*iter)->SetNotInSubgraph();
			(*iter)->SetDE( false );
			int pos = (*iter)->GetPositionOfContig( c );
			int ori = (*iter)->GetOrientationOfContig( c );
			//int preOri = ori;
			if( !isFirst )
				(*iter)->ReplaceContig( pos, m_contigsArray[ scaffold->GetID() ] );

			if( !ifSame ){
				if( ori == MINUS )
					ori = PLUS;
				else
					ori = MINUS;
				
				(*iter)->SetOri( pos, ori );
			}

			/*if( (*iter)->GetOriString().find( "1437\t" ) != string::npos ){
				cerr<<"find update position\n";
				cerr<<(*iter)->GetOriString()<<endl;
				cerr<<"c: "<<c->GetName()<<endl;
				cerr<<"ifSame: "<<ifSame<<endl;
				cerr<<"ori: ";
				if( ori == PLUS )
					cerr<<"PLUS\n";
				else
					cerr<<"MINUS\n";
				cerr<<"Pos: ";
				if( pos == START )
					cerr<<"START\n";
				else
					cerr<<"END\n";
				cerr<<"leftDis: "<<c->GetLeftDis()<<endl;
				cerr<<"rightDis: "<<c->GetRightDis()<<endl;
				cerr<<endl;
			}
			*/

			// recalculate the distance
#ifdef CHECK
			int oldDis = (*iter)->GetDis();
#endif

			// FIXME: if two contigs have the same repetitive repeat, add the length of the repeat to the distance
			if( pos == START ){
				if( ori == PLUS ){
						(*iter)->SetDis( (*iter)->GetDis() - (int) c->GetRightDis() );
				}
				else{
						(*iter)->SetDis( (*iter)->GetDis() - (int)c->GetLeftDis() );
				}
			}
			else{
				if( ori == PLUS ){
						(*iter)->SetDis( (*iter)->GetDis() - (int) c->GetLeftDis() );
				}
				else{
						(*iter)->SetDis( (*iter)->GetDis() - (int) c->GetRightDis() );
				}
			}

			//cout<<"after: "<<(*iter)->ToString()<<endl;

#ifdef CHECK
			if( (*iter)->GetDis() > oldDis ){
				//cout<<"new distance > old distance"<<endl;
				//cout<<"edge is "<<(*iter)->GetOriString()<<endl<<endl;
			}
#endif

			// check if this edge distance is happy or not                                                                                                                                                                                                        
                        // if not happy, delete this edge                                                                                                                                                                                                                     
                        if( (*iter)->GetDis() + (*iter)->GetStd() * Configure::STD_TIMES < -Configure::KMER ){
#ifdef DEBUG
                                cerr<<"Try to delete edges from other contig: "<<(*iter)->GetOriString()<<endl;
#endif
                                // remove from the other contig                                                                                                                                                                                                               
                                //Contig *tempOtherContig = (*iter)->GetOtherContig( c );                                                                                                                                                                                     
#ifdef DEBUG
                                cerr<<"other contig is "<<tempOtherContig->GetName()<<endl;
#endif

                                // locate the edge                                                                                                                                                                                                                            
                                bool findEdge = false;
                                list<PET*> *tempEdges = tempOtherContig->GetLeftEdges();
                                for( list<PET*>::iterator tempIter = tempEdges->begin(); tempIter != tempEdges->end(); tempIter++ ){
#ifdef DEBUG
                                        cerr<<"checking "<<(*tempIter)->GetOriString()<<endl;
#endif
                                        if( *tempIter == *iter ){
#ifdef DEBUG
                                                cerr<<"find"<<endl;
#endif
                                                tempEdges->erase( tempIter );
                                                findEdge = true;
                                                break;
                                        }
                                }
                                if( !findEdge ){
                                        tempEdges = tempOtherContig->GetRightEdges();
                                        for( list<PET*>::iterator tempIter = tempEdges->begin(); tempIter != tempEdges->end(); tempIter++ ){
#ifdef DEBUG
                                                cerr<<"checking "<<(*tempIter)->GetOriString()<<endl;
#endif
                                                if( *tempIter == *iter ){
#ifdef DEBUG
                                                        cerr<<"find"<<endl;
#endif
                                                        tempEdges->erase( tempIter );
                                                        findEdge = true;
                                                        break;
						}
                                        }

                                }

#ifdef DEBUG
                                cerr<<"deleted from other contig? "<<findEdge<<endl;
                                cerr<<"remove from current contig"<<endl;
#endif

                                // remove from current contig                                                                                                                                                                                                                 
                                PET *delEdge = *iter;
                                iter = edges->erase( iter );
                                delete delEdge;
                                continue;
#ifdef DEBUG
                                cerr<<"done"<<endl;
#endif
                        }


			if( !isFirst ){
				// add current edge to super contig
				m_contigsArray[ scaffold->GetID() ]->AddEdge( *iter );

				// remove current edge from current contig
				iter = edges->erase( iter );
			}
			else
				iter++;
		}
	}
}

// release the memory of edges
void Graph::DeleteEdgesMultiLib( Contig *c, list<PET*> *edges, bool isFirst, ScaffoldResult *scaffold, bool ifSame, FILE* discardEdgeFile ){
	bool isSpecial = false;
	

	list<PET*>::iterator iter = edges->begin();
	while( iter != edges->end() ){
		//cout<<"update edge\t"<<(*iter)->GetOriString()<<endl;
		if( isSpecial )
			cerr<<"New Edge\n";
		
		// if this edge is within the subgraph
		if( (*iter)->GetStartContig()->isInSubgraph() && (*iter)->GetEndContig()->isInSubgraph()
		    && ( (!(*iter)->GetStartContig()->IsRepeat() && !(*iter)->GetEndContig()->IsRepeat() 
			  && (*iter)->GetStartContig()->GetScaffoldID() == (*iter)->GetEndContig()->GetScaffoldID() ) 
			 || (*iter)->GetStartContig()->IsRepeat() || (*iter)->GetEndContig()->IsRepeat() ) ){
			
		    //&& (*iter)->GetStartContig()->GetScaffoldID() == (*iter)->GetEndContig()->GetScaffoldID() ){

			if( isSpecial )
				cerr<<"delete edges: "<<(*iter)->GetStartContig()->GetName()<<"\t"<<(*iter)->GetEndContig()->GetName()<<endl;
			
			if( !(*iter)->IfVisited() ){
				Contig *otherContig = (*iter)->GetOtherContig( c );
				int otherPos = (*iter)->GetPositionOfContig( otherContig );
				int otherOri = (*iter)->GetOrientationOfContig( otherContig );
				if( !otherContig->IsRepeat() ){
					// if it is the first time, update the status of edge
					(*iter)->VisitEdge();
					iter = edges->erase( iter );
				}
				else{
					// directly delete this edge
					PET *temp = *iter;
					iter = edges->erase( iter );
					// delete from repeat
					if( (otherOri == PLUS && otherPos == START) || (otherOri == MINUS && otherPos == END) ){
						// right edge
						otherContig->RemoveOriginalRepetitiveEdgeMultiLib( temp, RIGHT );
					}
					else{
						// left edge
						otherContig->RemoveOriginalRepetitiveEdgeMultiLib( temp, LEFT );
					}
					// remove from the PET library
					temp->RemoveFromMultiLib();
					delete temp;
				}
				// if it is the first time, update the status of edge
				//(*iter)->VisitEdge();
				//fprintf( discardEdgeFile, "%s\n", (*iter)->GetOriString().c_str() );
				//iter = edges->erase( iter );

			}
			else{
				if( isSpecial )
					cerr<<"second time\n";
				// the second time, delete
				PET *temp = *iter;
				//cout<<"delete: "<<temp->GetOriString()<<endl;
				// delete from library
				(*iter)->DeleteFromLib();

				iter = edges->erase( iter );
	
				delete temp;
				if( isSpecial )
					cerr<<"second time done\n";
			}
			if( isSpecial )
				cerr<<"delete edges succeed\n";
		}
		else if( c->GetScaffoldID() != -1){
			if( isSpecial )
				cerr<<"update edges\n";

			/*if( (*iter)->GetStartContig()->GetName() == "contig_42"  && (*iter)->GetEndContig()->GetName() == "contig_43" ){
				cout<<"find the edge\n";
				}*/

			//cout<<"scaffold:\n"<<scaffold->GetScaffoldString()<<endl;
			//cout<<"contig:\n"<<c->GetScaffoldString()<<endl;
			//cout<<"Edge is needed to be updated"<<endl;
			//cout<<"before: "<<(*iter)->ToString()<<endl;
			
			// if this edge connects outside the subgraph, update it
			/*if( isFirst ){
				// meet the first contig in a scaffold, replace it
				if( scaffold->GetID() == -1 ){
					scaffold->SetID( c->GetID() );
					// update length of contig
					c->SetLength( scaffold->GetLength() );
				}
			}*/

			// replace edges
			int pos = (*iter)->GetPositionOfContig( c );
			int ori = (*iter)->GetOrientationOfContig( c );
			//int preOri = ori;
			if( !isFirst )
				(*iter)->ReplaceContig( pos, m_contigsArray[ scaffold->GetID() ] );

			if( !ifSame ){
				if( ori == MINUS )
					ori = PLUS;
				else
					ori = MINUS;
				
				(*iter)->SetOri( pos, ori );
			}

			// recalculate the distance
#ifdef CHECK
			int oldDis = (*iter)->GetDis();
#endif

			// FIXME: if two contigs have the same repeat at the end, remove the length of this repeat from the distance
			if( pos == START ){
				if( ori == PLUS  ){
						(*iter)->SetDis( (*iter)->GetDis() - (int) c->GetRightDis() );
						(*iter)->SetDisWithGap( (*iter)->GetDisWithGap() - (int) c->GetRightDisWithGap() );
				}
				else{
						(*iter)->SetDis( (*iter)->GetDis() - (int)c->GetLeftDis() );
						(*iter)->SetDisWithGap( (*iter)->GetDisWithGap() - (int)c->GetLeftDisWithGap() );
				}
			}
			else{
				if( ori == PLUS ){
						(*iter)->SetDis( (*iter)->GetDis() - (int) c->GetLeftDis() );
						(*iter)->SetDisWithGap( (*iter)->GetDisWithGap() - (int) c->GetLeftDisWithGap() );
				}
				else{
						(*iter)->SetDis( (*iter)->GetDis() - (int) c->GetRightDis() );
						(*iter)->SetDisWithGap( (*iter)->GetDisWithGap() - (int) c->GetRightDisWithGap() );
				}
			}

			//cout<<"without gap: "<<(*iter)->GetDis()<<endl;
			//cout<<"with gap: "<<(*iter)->GetDisWithGap()<<endl;

			//cout<<"after: "<<(*iter)->ToString()<<endl<<endl;
			
#ifdef CHECK
			if( (*iter)->GetDisWithGap() > oldDis ){
				cout<<"new distance > old distance: "<<(*iter)->GetDisWithGap()<<"\t"<<oldDis<<endl;
				cout<<"edge is "<<(*iter)->GetOriString()<<endl<<endl;
			}
#endif

			if( !isFirst ){
				// add current edge to super contig
				m_contigsArray[ scaffold->GetID() ]->AddEdgeMultiLib( *iter );

				// remove current edge from current contig
				iter = edges->erase( iter );
			}
			else
				iter++;
			
			if( isSpecial )
				cerr<<"update edges succeed\n";
		}
	}

	if( isSpecial )
		cerr<<"finish updage edges of multiple libraries\n";
}

// check if current scaffold has the same orientation in scaffold and edge
bool Graph::IfSameOriInScaffoldAndEdge( Contig *c, int ori, ScaffoldResult *s ){
	/*string oriString = "BE";
	if( ori == MINUS )
		oriString = "EB";*/

	string contigName = c->GetNameOfFirstContig();
	string oriString = c->GetOriOfFirstContig( ori );

	vector<string> lines;
	Split( s->GetScaffoldString(), "\n", &lines );

	for( int i = 0; i < (int) lines.size(); i++ ){
		vector<string> line;
		Split( lines.at( i ), "\t", &line );
		if( line.at( 0 ) == contigName ){
			if( line.at( 1 ) == oriString )
				return true;
			else 
				return false;
		}
	}

	return true;
}


// find subgraph
list<Contig*> * Graph::FindSubgraph( int &numOfContig, int &numOfBorderContig, int &minClusterSize ){
	list<Contig*> *subgraph = NULL;
		
	Contig *c = NULL;

	bool hasBorderContig = true;

	// check if there is any border contigs
	//cerr<<"total border contigs: "<<m_borderContigsList->size()<<endl;
	if( !m_borderContigsList->empty() ){
		for( list<Contig*>::iterator iter = m_borderContigsList->begin(); iter != m_borderContigsList->end(); iter++ ){
			if( !(*iter)->IsRepeat() ){
				c = *iter;
				break;
			}
		}

		if( c == NULL )
			hasBorderContig = false;
		else{
			//cerr<<"Starting contig is: "<<c->GetName()<<endl;
			if( c->HasRightEdge() ){
				//cerr<<"traverse right\n";
				//subgraph = TraverseSubgraph( *(m_borderContigsList->begin()), RIGHT, 
				subgraph = TraverseSubgraph( c, RIGHT, 
							     numOfContig, numOfBorderContig, minClusterSize );
			}
			else{
				//cerr<<"traverse left\n";
				//subgraph = TraverseSubgraph( *(m_borderContigsList->begin()), LEFT, 
				subgraph = TraverseSubgraph( c, LEFT, 
							     numOfContig, numOfBorderContig, minClusterSize );
			}
		}
	}
	else
		hasBorderContig = false;
	
	if( !hasBorderContig && !m_contigsList->empty() ){
		for( list<Contig*>::iterator iter = m_contigsList->begin(); iter != m_contigsList->end(); iter++ ){
			if( !(*iter)->IsRepeat() ){
				c = *iter;
			}
		}
		//c = *( m_contigsList->begin() );
		//subgraph = TraverseSubgraph( *(m_contigsList->begin()), BOTH, 
		subgraph = TraverseSubgraph( c, BOTH, 
					     numOfContig, numOfBorderContig, minClusterSize );
		/*if( c->HasRightEdge() )
			subgraph = TraverseSubgraph( *(m_contigsList->begin()), RIGHT, 
						numOfContig, numOfBorderContig, minClusterSize );
		else
			subgraph = TraverseSubgraph( *(m_contigsList->begin()), LEFT, 
					numOfContig, numOfBorderContig, minClusterSize );
		*/
	}

	// check if need to add edges
	for( list<Contig*>::iterator iter = subgraph->begin(); iter != subgraph->end(); iter++ ){
		list<PET*> *edges = (*iter)->GetLeftEdges();
		for( list<PET*>::iterator edgeIter = edges->begin(); edgeIter != edges->end(); edgeIter++ ){
			// don't add new edges for repeats to avoid adding edges with wrong orientations
			if( (*edgeIter)->GetStartContig()->IsRepeat() || (*edgeIter)->GetEndContig()->IsRepeat() )
				continue;

			if( (*edgeIter)->GetStartContig()->isInSubgraph() && (*edgeIter)->GetEndContig()->isInSubgraph() ){
				if( !(*edgeIter)->isInSubgraph() ){
					(*edgeIter)->SetInSubgraph();
					if( (*edgeIter)->GetSize() < minClusterSize )
						minClusterSize = (*edgeIter)->GetSize();
				}
			}
		}

		edges = (*iter)->GetRightEdges();
		for( list<PET*>::iterator edgeIter = edges->begin(); edgeIter != edges->end(); edgeIter++ ){
			// don't add new edges for repeats to avoid adding edges with wrong orientations
			if( (*edgeIter)->GetStartContig()->IsRepeat() || (*edgeIter)->GetEndContig()->IsRepeat() )
				continue;

			if( (*edgeIter)->GetStartContig()->isInSubgraph() && (*edgeIter)->GetEndContig()->isInSubgraph() ){
				if( !(*edgeIter)->isInSubgraph() ){
					(*edgeIter)->SetInSubgraph();
					if( (*edgeIter)->GetSize() < minClusterSize )
						minClusterSize = (*edgeIter)->GetSize();
				}
			}
		}
	}

	return subgraph;
}

// count all the subgraph
void Graph::CountSubgraph(){
	list<Contig*> *subgraph = new list<Contig*>;
		
	Contig *c;
	int numOfContig;
	int numOfBorderContig;
	int numberOfRepeats;
	int minClusterSize;

	list<Contig*>::iterator iter = m_borderContigsList->begin();
	while( iter != m_borderContigsList->end() ){
		c = *iter;
		//cerr<<c->GetName()<<endl;
		
		subgraph->clear();
		numOfContig = 0;
		numOfBorderContig = 0;
		numberOfRepeats = 0;
		minClusterSize = 0;
		subgraph = TraverseSubgraphToCount( c, RIGHT, 
						    numOfContig, numOfBorderContig, minClusterSize, numberOfRepeats );
		//cerr<<"finish traversing to the right\n";
		AddEdgeToSubgraph( subgraph );
		//cerr<<"finish adding edges\n";
		//	if( subgraph->size() > 1 )
		//cerr<<"total contig in subgraph\t"<<numOfContig<<"\trepeats\t"<<numberOfRepeats<<endl;
		ClearSubgraph( subgraph );
		
		subgraph->clear();
		numOfContig = 0;
		numOfBorderContig = 0;
		numberOfRepeats = 0;
		minClusterSize = 0;
		subgraph = TraverseSubgraphToCount( c, LEFT, 
						    numOfContig, numOfBorderContig, minClusterSize, numberOfRepeats );
		//cerr<<"finish traversing to the left\n";
		AddEdgeToSubgraph( subgraph );
		//cerr<<"finish adding edges\n";
		//if( subgraph->size() > 1 )
		//	cerr<<"total contig in subgraph\t"<<numOfContig<<"\trepeats\t"<<numberOfRepeats<<endl;
		ClearSubgraph( subgraph );

		iter++;
		
	}
	
	iter = m_contigsList->begin();
	while( iter != m_contigsList->end() ){
		c = *iter;
		subgraph->clear();
		numOfContig = 0;
		numOfBorderContig = 0;
		numberOfRepeats = 0;
		minClusterSize = 0;
		subgraph = TraverseSubgraphToCount( c, BOTH, 
						    numOfContig, numOfBorderContig, minClusterSize, numberOfRepeats );
		AddEdgeToSubgraph( subgraph );
		//if( subgraph->size() > 1 )
		//	cerr<<"total contig in subgraph\t"<<numOfContig<<"\trepeats\t"<<numberOfRepeats<<endl;
		ClearSubgraph( subgraph );
	       

		iter++;
	}

	subgraph->clear();
	delete subgraph;

}

// unlabel all the contigs in the subgraph
void Graph::ClearSubgraph( list<Contig*> *subgraph ){
	list<Contig*>::iterator iter = subgraph->begin();
	while( iter != subgraph->end() ){
		(*iter)->SetNotInSubgraph();
		(*iter)->SetExtensionOri( NEITHER );
		iter++;
	}
}

// check if there are missing edges in this subgraph
void Graph::AddEdgeToSubgraph( list<Contig*> *subgraph ){
	for( list<Contig*>::iterator iter = subgraph->begin(); iter != subgraph->end(); iter++ ){
		list<PET*> *edges = (*iter)->GetLeftEdges();
		for( list<PET*>::iterator edgeIter = edges->begin(); edgeIter != edges->end(); edgeIter++ ){
			if( (*edgeIter)->GetStartContig()->isInSubgraph() && (*edgeIter)->GetEndContig()->isInSubgraph() ){
				if( !(*edgeIter)->isInSubgraph() ){
					(*edgeIter)->SetInSubgraph();
				}
			}
		}

		edges = (*iter)->GetRightEdges();
		for( list<PET*>::iterator edgeIter = edges->begin(); edgeIter != edges->end(); edgeIter++ ){
			if( (*edgeIter)->GetStartContig()->isInSubgraph() && (*edgeIter)->GetEndContig()->isInSubgraph() ){
				if( !(*edgeIter)->isInSubgraph() ){
					(*edgeIter)->SetInSubgraph();
				}
			}
		}
	}
}

// traverse and find subgraph starting from contig c 
// c is the start contig; direction is the extension direction
list<Contig*> * Graph::TraverseSubgraphToCount( Contig *c, int direction, int &numOfContig, int &numOfBorderContig, int &minClusterSize, int &numberOfRepeats ){
	list<Contig*> *subgraph = new list<Contig*>;

	numOfContig = 0;
	numOfBorderContig = 0;
		
	list<Contig*> possibleContigs;
	possibleContigs.push_back( c );
	c->SetExtensionOri( direction );

	while( !possibleContigs.empty() ){
		// get the first contig
		Contig *currentContig = *possibleContigs.begin();
		possibleContigs.pop_front();

		// label this contig as in subgraph if it is not
		if( !currentContig->isInSubgraph() ){
			currentContig->SetInSubgraph();
			subgraph->push_back( currentContig );
			//cerr<<"subgraph: "<<currentContig->GetName()<<"\t"<<currentContig->GetContigType()<<endl;
			numOfContig++;
			if( currentContig->IsBorderContig() )
				numOfBorderContig++;
			if( currentContig->GetContigType() == REPEAT )
				numberOfRepeats++;
		}

		// traverse edges
		if( currentContig->GetExtensionOri() == RIGHT )
			AddEdgesToSubgraph( currentContig, currentContig->GetRightEdges(), &possibleContigs, minClusterSize );
		else if( currentContig->GetExtensionOri() == LEFT )
			AddEdgesToSubgraph( currentContig, currentContig->GetLeftEdges(), &possibleContigs, minClusterSize );
		else if( currentContig->GetExtensionOri() == BOTH ){
			AddEdgesToSubgraph( currentContig, currentContig->GetLeftEdges(), &possibleContigs, minClusterSize );
			AddEdgesToSubgraph( currentContig, currentContig->GetRightEdges(), &possibleContigs, minClusterSize );
		}
	}

	return subgraph;
}

int Graph::GetNumOfContigs(){
	return m_numberOfContigs;
}

// traverse and find subgraph starting from contig c 
// c is the start contig; direction is the extension direction
list<Contig*> * Graph::TraverseSubgraph( Contig *c, int direction, int &numOfContig, int &numOfBorderContig, int &minClusterSize ){
	list<Contig*> *subgraph = new list<Contig*>;

	numOfContig = 0;
	numOfBorderContig = 0;
		
	list<Contig*> possibleContigs;
	possibleContigs.push_back( c );
	c->SetExtensionOri( direction );
	
	//cerr<<"Find subgraph\n";

	while( !possibleContigs.empty() ){
		// get the first contig
		Contig *currentContig = *possibleContigs.begin();
		possibleContigs.pop_front();

		// label this contig as in subgraph if it is not
		if( !currentContig->isInSubgraph() ){
			//cerr<<"add contig to subgraph:\t"<<currentContig->GetName()<<endl;
			currentContig->SetInSubgraph();
			subgraph->push_back( currentContig );
			numOfContig++;
			if( currentContig->IsBorderContig() )
				numOfBorderContig++;
		}

		if( currentContig->IsRepeat() ){
			continue;
		}

		// traverse edges
		if( currentContig->GetExtensionOri() == RIGHT )
			AddEdgesToSubgraph( currentContig, currentContig->GetRightEdges(), &possibleContigs, minClusterSize );
		else if( currentContig->GetExtensionOri() == LEFT )
			AddEdgesToSubgraph( currentContig, currentContig->GetLeftEdges(), &possibleContigs, minClusterSize );
		else if( currentContig->GetExtensionOri() == BOTH ){
			AddEdgesToSubgraph( currentContig, currentContig->GetLeftEdges(), &possibleContigs, minClusterSize );
			AddEdgesToSubgraph( currentContig, currentContig->GetRightEdges(), &possibleContigs, minClusterSize );
		}
	}

	//cerr<<endl;
	return subgraph;
}


// add list of edges into subgraph
// add edges to subgraph, and all contigs except c to contigList
void Graph::AddEdgesToSubgraph( Contig *c, list<PET*> *edges, list<Contig*> *contigList, int &minClusterSize ){
	list<PET*>::iterator iter;
	for( iter = edges->begin(); iter != edges->end(); iter++ ){
		if( (*iter)->isInSubgraph() )
			continue;
		else{
			if( (*iter)->GetSize() < minClusterSize )
				minClusterSize = (*iter)->GetSize();

			//cerr<<"Add edge: "<<(*iter)->ToString()<<endl;
			(*iter)->SetInSubgraph();
			(*iter)->SetDE( false );
			Contig *otherContig = (*iter)->GetOtherContig( c );
			int otherPos = (*iter)->GetPositionOfContig( otherContig );
			int otherOri = (*iter)->GetOrientationOfContig( otherContig );
			//int ext;
			if( otherContig->GetLength() < Configure::UPPERBOUND || otherContig->GetContigType() == REPEAT ){
				// short contig or repeats, just add and extend to both orientation
				otherContig->SetExtensionOri( BOTH );
			}
			else{
				// border contig, need consider orientation
				if( ( otherPos == START && otherOri == MINUS ) 
				|| ( otherPos == END && otherOri == PLUS ) ){
					// new contig should extend to left
					otherContig->SetExtensionOri( LEFT );
				}
				else
					otherContig->SetExtensionOri( RIGHT );
			}

			contigList->push_back( otherContig );
		}
	}
}

// check if current graph still has edges
bool Graph::HasEdges(){
	return !m_contigsList->empty();
}

// output final scaffolds to files
// return -1 if failed
int Graph::OutputScaffolds( string fileName ){
	ofstream scafWriter( fileName.c_str() );
	if( scafWriter == NULL ){
		cout<<"ERROR: Cannot open "<<fileName<<" file"<<endl;
		return -1;
	}

	// sort scaffold
	m_scaffoldsList->sort( compare_scaffold );
	//cerr<<"finish sorting scaffold\n";

	int num = 1;
	for( list<Contig*>::iterator iter = m_scaffoldsList->begin(); iter != m_scaffoldsList->end(); iter++ ){
		string scaf = (*iter)->GenScaffoldString( PLUS ) + "0\n";
		char buffer[ 200 ];
		sprintf( buffer, ">%s%d\tlength: %.0f\tcov: %.1f\n", Configure::SCAFFOLD_PREFEX.c_str(), num, (*iter)->GetLength(), (*iter)->GetCov() );
		string head( buffer );
		scafWriter.write( head.c_str(), head.length() );
		scafWriter.write( scaf.c_str(), scaf.length() );
		num++;
	}
	scafWriter.close();
	//cerr<<"finish wrting scaffolds\n";
	return 1;
}

// comparison, not case sensitive.
bool compare_scaffold( Contig *first, Contig *second)
{
	if( first->GetLength() > second->GetLength() )
		return true;
	else 
		return false;
}

// generate the id map of all contigs
void Graph::GenerateIDMap()
{
	if( Configure::FILE_TYPE == VELVET ){
		// handle velvet format
		for( int i = 0; i < m_numberOfContigs; i++ ){
			vector<string> *name = new vector<string>;
			Split( m_contigsArray[ i ]->GetName(), "_", name );
			m_contigID->insert( pair<int, Contig*>( atoi( name->at( 1 ).c_str() ), m_contigsArray[ i ] ) );
		}
	}
}

// get the contig using velvet id
Contig* Graph::GetContigUsingID( int id )
{
	return (*m_contigID)[ id ];
}

// clear the edges number
void Graph::ClearEdge(){
	m_numberOfEdges = 0;
}

// remove one edge from graph
void Graph::RemoveOneEdge()
{
	this->m_numberOfEdges--;
}

// set a contig as singleton
void Graph::SetAsSingleton( Contig *c )
{
	if( !c->HasEdge() )
	{
		// remove to singleton list
		m_scaffoldsList->push_back( c );
		m_contigsArray[ c->GetID() ] = NULL;
		if( c->IsBorderContig() )
			m_borderContigsList->erase( c->GetBorderListPos() );
		m_contigsList->erase( c->GetListPos() );
	}
}

// get all contigs
list<Contig*>* Graph::GetContigs()
{
	return m_contigsList;
}

// insert repeat contig name
void Graph::InsertRepeatContigName( Contig *c ){
	m_repeatContigSet->insert( pair<const char*, bool> (c->GetName().c_str(), true ) );
}

// check if a certain contig is a repeat
bool Graph::IsRepeat( string name ){
	hash_map<const char*, bool, hash<const char*>, eqName>::iterator pos = m_repeatContigSet->find( name.c_str() );
	if( pos == m_repeatContigSet->end() )
		return false;
	else
		return true;
}


// check the number of edges of certain contig
void Graph::CheckNumberOfEdges( string name ){
	for( list<Contig*>::iterator iter = m_borderContigsList->begin(); iter != m_borderContigsList->end(); iter++ ){
		if( (*iter)->GetName() == name ){
			cerr<<"find "<<name<<endl;
			cerr<<"scaffolds: "<<endl<<(*iter)->GetScaffoldString()<<endl;
			cerr<<"left edges: "<<(*iter)->GetLeftEdges()->size()<<endl;
			for( list<PET*>::iterator edgeIter = (*iter)->GetLeftEdges()->begin(); edgeIter != (*iter)->GetLeftEdges()->end(); edgeIter++ ){
				cerr<<"\t"<<(*edgeIter)->GetOriString()<<endl;
			}
			cerr<<"right edges: "<<(*iter)->GetRightEdges()->size()<<endl;
			for( list<PET*>::iterator edgeIter = (*iter)->GetRightEdges()->begin(); edgeIter != (*iter)->GetRightEdges()->end(); edgeIter++ ){
				cerr<<"\t"<<(*edgeIter)->GetOriString()<<endl;
			}
			cerr<<endl<<endl;
		}
	}
}
