#include "ContigConverter.h"

ContigConverter::ContigConverter(void)
{
	myContigs = new vector<Contig*>;
	m_repeatContigs = new list<Contig*>;
	m_smallContigs = new list<Contig*>;

	m_contigNameHashMap = new hash_map<const char*, int, hash<const char*>, eqName>;
	
	m_libString = "";
}

ContigConverter::~ContigConverter(void)
{
	myContigs->clear();
	delete myContigs;
	//DeleteContigs( m_repeatContigs );
	m_repeatContigs->clear();
	delete m_repeatContigs;
	//DeleteContigs( m_smallContigs );
	m_smallContigs->clear();
	delete m_smallContigs;
	m_contigNameHashMap->clear();
	delete m_contigNameHashMap;
}

void ContigConverter::DeleteContigs( list<Contig*> *contigs ){
	list<Contig*>::iterator iter = contigs->begin();
	while( iter != contigs->end() ){
		Contig* temp = *iter;
		iter = contigs->erase( iter );
		delete temp;
	}
}

// Convert contig file
// return: -1 if failed
int ContigConverter::ConvertContigFile( string fileName, Graph *graph, list<PetLibrary*> *libs ){
	ifstream contigReader( fileName.c_str() );

	if( contigReader == NULL )
	{
		cout<<"ERROR: Cannot open "<<fileName<<" file"<<endl;
		return -1;
	}

	if( Configure::FILE_FORMAT == FASTA ){
		if( Configure::FILE_TYPE == VELVET ){
			// analyze velvet file
			if( AnalyzeVelvet( &contigReader, myContigs ) == -1 )
				return -1;
		}
		else if( Configure::FILE_TYPE == SOAP ){
			// analyze soapdenovo file
			if( AnalyzeSoapDenovo( &contigReader, myContigs ) == -1 )
				return -1;
		}
		else if( Configure::FILE_TYPE == SOAP_CONTIG ){
			// analyze soapdenovo contig format
			//cerr<<"soap contig format\n";
			if( AnalyzeSoapDenovoContig( &contigReader, myContigs ) == -1 )
				return -1;
		}
		else if( Configure::FILE_TYPE == OPERA_SCAF ){
			//cerr<<"opera scaffold format\n";
			if( AnalyzeOperaScaf( &contigReader, myContigs ) == -1 )
				return -1;
		}
		else if( Configure::FILE_TYPE == NORMAL ){
			// analyze normal fasta file
			//cerr<<"normal type\n";
			if( AnalyzeNormalFasta( &contigReader, myContigs ) == -1 )
				return -1;
		}
	}
	else if( Configure::FILE_FORMAT == STATISTIC ){
		// read the contig file directly in opera format
		if( AnalyzeStatistic( &contigReader, myContigs ) == -1 )
			return -1;
	}

	// save all contig names in graph
	for( int i = 0; i < (int) myContigs->size(); i++ ){
		graph->SaveContigName( myContigs->at( i )->GetName() );
	}

	//cerr<<"Done saving all contig names\n";

	// calculate coverage using mapping information
	//if( Configure::FILE_TYPE == NORMAL && Configure::FILE_FORMAT != STATISTIC ){
	//if( Configure::MAP_TYPE != OPERA && Configure::FILE_FORMAT != STATISTIC ){
	if( Configure::FILE_FORMAT != STATISTIC ){
		if( CalculateCovUsingMapping( libs ) == -1 )
			return -1;
	}
	//}

	// update contig threshold, only use pair-end libraries
	if( 2 * Configure::MIN_LIBRARY_MEAN > 500 && 2 * Configure::MIN_LIBRARY_MEAN < 2000 ){
		Configure::CONTIG_SIZE_THERSHOLD = 2 * Configure::MIN_LIBRARY_MEAN;
	}

	cerr<<"min library mean : "<<Configure::MIN_LIBRARY_MEAN<<endl;
	cerr<<"minimum contig length is "<<Configure::CONTIG_SIZE_THERSHOLD<<endl;

	// if need to remove repeat, remove them
	//if( Configure::FILTER_REPEAT ){
	FilterRepeat( myContigs, graph );
	/*}
	  else{
	  // save the contigs directly and remove the small contigs
	  RemoveSmallContigs( myContigs, graph );
	  }
	*/

	// generate the ID map of all contigs
	graph->GenerateIDMap();

	// print repeat and small contigs
	if( PrintContigs( m_repeatContigs, Configure::OUTPUT_FOLDER + "repeatContigs" ) == -1 )
		return -1;

	if( PrintContigs( m_smallContigs, Configure::OUTPUT_FOLDER + "smallContigs" ) == -1 )
		return -1;

	// output the opera format contig file
	//if( Configure::FILE_FORMAT != STATISTIC ){
	if( graph->OutputContigs( Configure::OUTPUT_FOLDER + "contigs" ) == -1 )
		return -1;
	//}

	return 1;
}

// calculate coverage using mapping information
int ContigConverter::CalculateCovUsingMapping( list<PetLibrary*> *libs ){
	// get all mapping files and handle them one by one
	for( int libID = 0; libID < (int) Configure::MULTI_LIB_INFO->size(); libID++ ){
		if(Configure::MULTI_LIB_INFO->at( libID )-> GetMapType() == SAM || Configure::MULTI_LIB_INFO->at( libID )-> GetMapType() == BOWTIE){
			// analyze all reads
			cerr<<"Analyzing "<<(libID+1)<<" library: ";
			string tempFileName = Configure::MULTI_LIB_INFO->at( libID )->GetFileName();
			cerr<<tempFileName<<endl;
			
			PetLibrary *newLib = new PetLibrary( Configure::MULTI_LIB_INFO->at( libID )->GetFileName() );
			libs->push_back( newLib );
			//If there is at leat a mapping file the contig coverage module can be runned
			if(Configure::MULTI_LIB_INFO->at( libID )->GetMapType() != OPERA)
				Configure::MAP_TYPE = Configure::MULTI_LIB_INFO->at( libID )->GetMapType();
			
			//cerr<<"MAP_TYPE"<<(int)Configure::MAP_TYPE<<endl;cin.get();

			// set library mean
			if( Configure::MULTI_LIB_INFO->at( libID )->GetMean() != -1 && Configure::MULTI_LIB_INFO->at( libID )->GetStd() != -1 )
				{
					newLib->SetMean( Configure::MULTI_LIB_INFO->at( libID )->GetMean() );
					newLib->SetStd( Configure::MULTI_LIB_INFO->at( libID )->GetStd() );
					newLib->SetOri( Configure::MULTI_LIB_INFO->at( libID )->GetOri() );
					Configure::READ_ORI = Configure::MULTI_LIB_INFO->at( libID )->GetOri();
					Configure::CALCULATE_LIB = false;
					Configure::CALCULATE_ORI = false;
				}
			else
				Configure::CALCULATE_LIB = true;
		
			if( Configure::MULTI_LIB_INFO->at( libID )->GetFileName().substr( tempFileName.length() - 4, 4 ) == ".bam" ){
				// it is a bam format
				tempFileName = Configure::OUTPUT_FOLDER + "temporary.sam";
				string cmd = "mkfifo " + tempFileName;
				if( system( cmd.c_str() ) );
				cmd  = Configure::SAMDIR.c_str() + std::string("samtools view ") + Configure::MULTI_LIB_INFO->at( libID )->GetFileName() + " > " + tempFileName + "&"; 
				if( system( cmd.c_str() ) );
			}

			ifstream mapReader( tempFileName.c_str() );

			if( mapReader == NULL ){
				cout<<"ERROR: Cannot open "<<Configure::MULTI_LIB_INFO->at( libID )->GetFileName()<<" file"<<endl;
			
				if( newLib->GetFileName().substr( newLib->GetFileName().length() - 4, 4 ) == ".bam" ){
					string cmd = "rm " + Configure::OUTPUT_FOLDER + "temporary.sam";
					if( system( cmd.c_str() ) );
				}

				return -1;
			}

			///////
			// variables for reading mapping files
			///////

			int oriNum[ 5 ];
			for( int i = 0; i < 5; i++ )
				oriNum[ i ] = 0;

			int numOfPet = 0;

			vector<double> distance;
		
			int num[ 5 ];
			int sum[ 5 ];
			vector<double> allDis[ 5 ];
			for( int i = 0; i < 5; i++ ){
				num[ i ] = 0;
				sum[ i ] = 0;
				allDis[ i ].resize( 0 );
			}

			vector<string> *preColumn = new vector<string>;
			vector<string> *nextColumn = new vector<string>;
			vector<string> *preNames = new vector<string>;
			vector<string> *nextNames = new vector<string>;

			char preName[ 200 ], nextName[ 200 ];
		
			double sumReadLength = 0;
			double numReadLength = 0;
		
			double tempAverage = 0;
			double tempNum = 0;
			bool useThreshold = false;
			//double contigThreshold = 0;

			string preLine, nextLine;
			while( getline( mapReader, preLine ) ){
				if( preLine.substr( 0, 1 ) != "@" )
					break;
			}
		
			while( getline( mapReader, nextLine ) ){
				try{     
					if( nextLine == "" )
						continue;

					sscanf( preLine.c_str(), "%s", preName );
					sscanf( nextLine.c_str(), "%s", nextName );
					string preNameString = string( preName );
					string nextNameString = string( nextName );
		
					Split( preLine, "\t", preColumn );
					Split( nextLine, "\t", nextColumn );
			
					Contig *firstContig = NULL;
					Contig *secondContig = NULL;
					int firstContigPos = 0;
					int secondContigPos = 0;
					string firstContigOri;
					string secondContigOri;
					int firstReadLength = 0;
					int secondReadLength = 0;

					if( IsPair( preNameString, nextNameString ) ){
						// is a pair and both reads are mapped
				
						// use both reads for coverage calculation
						if( !(atoi(preColumn->at( 1 ).c_str()) & 0X4) ){
							// if the first read is mappable
							string contigName = preColumn->at( 2 );
				
							hash_map<const char*, int, hash<const char*>, eqName>::iterator pos = m_contigNameHashMap->find( contigName.c_str() );
							if( pos != m_contigNameHashMap->end() ){
								// add 1 to the number of reads in this contig
								Contig *currentContig = myContigs->at(  pos->second );
								int readLength = 0;
								if( Configure::MULTI_LIB_INFO->at( libID )->GetMapType() == BOWTIE )
									readLength = preColumn->at( 4 ).length();
								else if( Configure::MULTI_LIB_INFO->at( libID )->GetMapType() == SAM ){
									readLength = preColumn->at( 9 ).length();
								}
						
								currentContig->AddOneMappedBases( readLength );

								firstContig = currentContig;
								firstContigPos = atoi( (*preColumn)[3].c_str() );
								firstContigOri = atoi( (*preColumn)[1].c_str() );
								firstReadLength = readLength;
							}
						}

						if( !(atoi( nextColumn->at( 1 ).c_str()) & 0X4) ){
							// if the second read is mappable
							string contigName = nextColumn->at( 2 );
							hash_map<const char*, int, hash<const char*>, eqName>::iterator pos = m_contigNameHashMap->find( contigName.c_str() );
							if( pos != m_contigNameHashMap->end() ){
								// add 1 to the number of reads in this contig
								Contig *currentContig = myContigs->at(  pos->second );
								int readLength = 0;
								if( Configure::MULTI_LIB_INFO->at( libID )->GetMapType() == BOWTIE )
									readLength = preColumn->at( 4 ).length();
								else if( Configure::MULTI_LIB_INFO->at( libID )->GetMapType() == SAM ){
									readLength = preColumn->at( 9 ).length();
								}
						
								currentContig->AddOneMappedBases( readLength );

								secondContig = currentContig;
								secondContigPos = atoi( (*nextColumn)[3].c_str() );
								secondContigOri = atoi( (*nextColumn)[1].c_str() );
								secondReadLength = readLength;
							}
						}

						// check if both reads are mappable
						if( (atoi(preColumn->at( 1 ).c_str()) & 0X4) || (atoi( nextColumn->at( 1 ).c_str()) & 0X4) )
							{ // one of the reads is not mappable	
								//cout<<"find unmappable read\n"<<preLine<<endl<<nextLine<<endl;
								while( getline( mapReader, preLine ) ){
									if( preLine != "" )
										break;
								}
								continue;
							}

						if( Configure::MAP_TYPE == SAM )
							{ // convert the sam format into bowtie format
								//cout<<"previous line: \n"<<preColumn->at( 1 )<<"\t"<<preColumn->at( 4 )<<endl;
								ConvertSamToBowtie( preColumn );
								//cout<<"new line: \n"<<preColumn->at( 1 )<<"\t"<<preColumn->at( 4 )<<endl;
					
								//cout<<"previous line: \n"<<nextColumn->at( 1 )<<"\t"<<nextColumn->at( 4 )<<endl;
								ConvertSamToBowtie( nextColumn );
								//cout<<"new line: \n"<<nextColumn->at( 1 )<<"\t"<<nextColumn->at( 4 )<<endl;

								firstContigOri = (*preColumn)[1];
								secondContigOri = (*nextColumn)[1].c_str();
							}

						sumReadLength += (*preColumn)[ 4 ].length();
						numReadLength++;
						sumReadLength += (*nextColumn)[ 4 ].length();
						numReadLength++;

						// for the name, only select the word before the first space
						(*preColumn)[ 0 ] = preNameString;
						(*nextColumn)[ 0 ] = nextNameString;
			
						if( (*preColumn)[ 2 ] == (*nextColumn)[ 2 ] )
							{
								// mapped to the same contig
								double dis = CalculateReadDisOnSameContig( preColumn, nextColumn );

								if( useThreshold ){
									// check orientation
									string firstOri = (*preColumn)[ 1 ];
									string secondOri = (*nextColumn)[ 1 ];
									int firstPos = atoi( (*preColumn)[ 3 ].c_str() );
									int secondPos = atoi( (*nextColumn)[ 3 ].c_str() );
						
									if( ( firstOri == "+" && secondOri == "+" && firstPos > secondPos ) || ( firstOri == "-" && secondOri == "-" && firstPos < secondPos )  ){
										// forward2: 2 ->...-> 1
										oriNum[ FORWARD2 ]++;
										sum[ FORWARD2 ] += dis;
										num[ FORWARD2 ]++;
										allDis[ FORWARD2 ].push_back( dis );
									}
									else if( ( firstOri == "+" && secondOri == "-" && firstPos < secondPos ) || ( firstOri == "-" && secondOri == "+" && firstPos > secondPos )  ){
										// in: -> <-
										oriNum[ IN ]++;
										sum[ IN ] += dis;
										num[ IN ]++;
										allDis[ IN ].push_back( dis );
									}
									else if( ( firstOri == "-" && secondOri == "+" && firstPos < secondPos ) || ( firstOri == "+" && secondOri == "-" && firstPos > secondPos )  ){
										// out: <- ->
										oriNum[ OUT ]++;
										sum[ OUT ] += dis;
										num[ OUT ]++;
										allDis[ OUT ].push_back( dis );
									}
									else if( ( firstOri == "+" && secondOri == "+" && firstPos < secondPos ) || ( firstOri == "-" && secondOri == "-" && firstPos > secondPos ) ){
										// forward1: 1 ->...-> 2
										oriNum[ FORWARD1 ]++;
										sum[ FORWARD1 ] += dis;
										num[ FORWARD1 ]++;
										allDis[ FORWARD1 ].push_back( dis );
									}
									else{
										oriNum[ ERROR ]++;
										sum[ ERROR ] += dis;
										num[ ERROR ]++;
										allDis[ ERROR ].push_back( dis );
									}
								}
					
								// check contig length
								hash_map<const char*, int, hash<const char*>, eqName>::iterator pos = m_contigNameHashMap->find( (*preColumn)[2].c_str() );
								Contig *c = NULL;
								if( pos != m_contigNameHashMap->end() )
									c = myContigs->at(  pos->second );

								if( !useThreshold || ( useThreshold && c != NULL && c->GetLength() > tempAverage * 3 ) )
									{
						
					
										// check the temporary average distance and contig length
										if( !useThreshold ){
											tempAverage = ( tempAverage * tempNum + dis ) / ( tempNum + 1 ); 
											tempNum++;
							
											if( tempNum >= 1000 ){
												useThreshold = true;
												//cerr<<"temp average is: "<<tempAverage<<endl;
											}
										}
									}
							}
						else{
							// map on different contigs, save this pair
							if( firstContig != NULL && secondContig != NULL ){
								newLib->AddMapInfo( firstContig, firstContigPos, firstContigOri, firstReadLength );
								newLib->AddMapInfo( secondContig, secondContigPos, secondContigOri, secondReadLength );
								numOfPet++;
							}
						}

						// get a new line
						while( getline( mapReader, preLine ) ){
							if( preLine != "" )
								break;
						}
					}
					else{
						// not a pair or at least one read is not mappable
				
						// use preLine to calculate coverage
						if( !(atoi(preColumn->at( 1 ).c_str()) & 0X4) ){
							// if the first read is mappable
							string contigName = preColumn->at( 2 );
				
							hash_map<const char*, int, hash<const char*>, eqName>::iterator pos = m_contigNameHashMap->find( contigName.c_str() );
							if( pos != m_contigNameHashMap->end() ){
								// add 1 to the number of reads in this contig
								Contig *currentContig = myContigs->at(  pos->second );
								int readLength = 0;
								if( Configure::MULTI_LIB_INFO->at( libID )->GetMapType() == BOWTIE )
									readLength = preColumn->at( 4 ).length();
								else if( Configure::MULTI_LIB_INFO->at( libID )->GetMapType() == SAM ){
									readLength = preColumn->at( 9 ).length();
								}
						
								currentContig->AddOneMappedBases( readLength );
							}
						}
				
						// move forward
						preLine = nextLine;
					}
				}
				catch( const std::out_of_range& oor) {
					std::cerr << "Out of Range error: " << oor.what() << '\n';
					std::cerr<<"preline is: "<<preLine<<endl;
					cerr<<"nextLine is: "<<nextLine<<endl;
				}
			}
		
			if( Configure::MULTI_LIB_INFO->at( libID )->GetFileName().substr( Configure::MULTI_LIB_INFO->at( libID )->GetFileName().length() - 4, 4 ) == ".bam" ){
				// it is a bam formattempFileName = Configure::OUTPUT_FOLDER + "temporary.sam";
				string cmd = "rm " + tempFileName;
				if( system( cmd.c_str() ) );
			}


			// calculate the library mean and std
			double avgReadLength = sumReadLength / numReadLength;
			newLib->SetReadLength( avgReadLength * 2 );
		
			int maxID = -1;
			if( Configure::CALCULATE_ORI ){
				int max = 0;
				maxID = -1;
				for( int i = 0; i < 5; i++ ){
					if( oriNum[ i ] > max ){
						max = oriNum[ i ];
						maxID = i;
					}
				}
			
				if( maxID == ERROR ){
					cout<<"ERROR: The orientation of reads is not \"in\", \"out\", neither \"forward\", Opera could not handle such orientation. "
					    <<"Please check the mapping files. Please feel free to contact gaosong@nus.edu.sg for further support. \n";
				
					if( newLib->GetFileName().substr( newLib->GetFileName().length() - 4, 4 ) == ".bam" ){
						string cmd = "rm " + Configure::OUTPUT_FOLDER + "temporary.sam";
						if( system( cmd.c_str() ) );
					}
				
					return -1;
				}
			
				newLib->SetOri( maxID );
				Configure::READ_ORI = maxID;
			}
		
			maxID = Configure::READ_ORI;
			//cerr<<"max id is "<<maxID<<endl;

			// calculate mean and std using IQR to exclude outliers
			double mean = 0;
			double std = 0;
			double percentageOfOutliers = 0;
			double numberOfUsedPairs = 0;

			if( Configure::CALCULATE_LIB ){
				vector<string> *path = new vector<string>;
				Split( newLib->GetFileName(), "/", path );
				vector<string> *names = new vector<string>;
				Split( path->at( path->size() - 1 ), ".", names );
				if( maxID != -1 ){
					if( IQR( &allDis[ maxID ], mean, std, names->at( 0 ), percentageOfOutliers, numberOfUsedPairs, newLib ) == -1 ){
						if( newLib->GetFileName().substr( newLib->GetFileName().length() - 4, 4 ) == ".bam" ){
							string cmd = "rm " + Configure::OUTPUT_FOLDER + "temporary.sam";
							if( system( cmd.c_str() ) );
						}
						path->clear();
						delete path;
						names->clear();
						delete names;
						return -1;
					}
				}
				path->clear();
				delete path;
				names->clear();
				delete names;

				newLib->SetMean( (int) mean );
				newLib->SetStd( (int) std );
				//cout<<"mean: "<<mean<<endl;
				//cout<<"std: "<<std<<endl;
			}
			else{
				// don't need to calculate the mean
				vector<string> *path = new vector<string>;
				Split( newLib->GetFileName(), "/", path );
				vector<string> *names = new vector<string>;
				Split( path->at( path->size() - 1 ), ".", names );

				ofstream disWriter( (Configure::OUTPUT_FOLDER + "distanceOfMatePairs_" + names->at( 0 )).c_str(), ios::out );
				if( disWriter == NULL ){
					cout<<"ERROR: Cannot open "<<(Configure::OUTPUT_FOLDER + "distanceOfMatePairs_" + names->at(0))<<" file"<<endl;
					path->clear();
					delete path;
					names->clear();
					delete names;
					return -1;
				}

				// get the threshold for outliers
				double lowerbound = newLib->GetMean() - Configure::STD_TIMES * newLib->GetStd();
				double upperbound = newLib->GetMean() + Configure::STD_TIMES * newLib->GetStd();
			
				// initialize the possibility of distance
				newLib->InitPossibility( lowerbound, upperbound );
			
				// calculate mean and std excluding all outliers
				double sum = 0;
				double num = 0;
				string disStr = "";
			
				for( int i = 0; i < (int) allDis[ maxID ].size(); i++ ){
					//cerr<<"i\n";
					disStr += itos( allDis[ maxID ].at( i ) ) + "\n";
					if( allDis[ maxID ].at( i ) >= lowerbound && allDis[ maxID ].at( i ) <= upperbound ){
						sum += allDis[ maxID ].at( i );
						newLib->AddOccuOfDistance( allDis[ maxID ].at( i ) );
						//distributionResults += itos( dis->at( i ) ) + "\n";
						num++;
					}
				}
			
				//cerr<<"start calculating pos\n";
				newLib->CalculatePossibility();
				//cerr<<"finish calculating pos\n";
			
				disWriter.write( disStr.c_str(), disStr.length() );
				percentageOfOutliers = (1.0 - num/(double)allDis[ maxID ].size()) * 100;
				numberOfUsedPairs = num;

				path->clear();
				delete path;
				names->clear();
				delete names;
				//cerr<<"finish writing distance\n";
			}
		
			// record the info of this library
			m_libString += "For the mapping file: " + newLib->GetFileName() + "\n";
			m_libString += "Mean length of the library is: " + itos( (int)newLib->GetMean() ) + "\n";
			m_libString.append( "Standard deviation of the library is: " + itos( (int)newLib->GetStd() ) + "\n" );
			m_libString.append( "Orientation of paired-reads is: " );
		
			if( Configure::READ_ORI == IN )
				m_libString.append( "in ->...<-\n" );
			else if( Configure::READ_ORI == OUT )
				m_libString.append( "out <-...->\n" );
			else if( Configure::READ_ORI == FORWARD2 )
				m_libString.append( "forward ->...->  2nd read...1st read\n" );
			else if( Configure::READ_ORI == FORWARD1 )
				m_libString.append( "forward ->...->  1st read...2nd read\n" );
		
			m_libString.append( "Number of used paired reads for estimating the mean and the standard deviation of the library is: "
					    + itos( numberOfUsedPairs ) + "\n" );
			char buffer[50];
			//int n, a=5, b=3;
			//n=sprintf (buffer, "Percentage of outliers is: %.2f\%\n", percentageOfOutliers );
			sprintf (buffer, "Percentage of outliers is: %.2f%%\n", percentageOfOutliers );
			string tempStr( buffer );
			m_libString.append( tempStr );
			m_libString.append( "\n" );

			// save minimum library mean
			if( newLib->GetMean() < Configure::MIN_LIBRARY_MEAN ){
				Configure::MIN_LIBRARY_MEAN = newLib->GetMean();
			}
		
			mapReader.close();
		
			preColumn->clear();
			delete preColumn;
			nextColumn->clear();
			delete nextColumn;
			preNames->clear();
			delete preNames;
			nextNames->clear();
			delete nextNames;
		
		}
	}

	// output library information
	ofstream libWriter( (Configure::OUTPUT_FOLDER + "lib.txt").c_str() );
	if( libWriter == NULL ){
		cout<<"ERROR: Cannot open "<<(Configure::OUTPUT_FOLDER + "lib.txt")<<" file"<<endl;
		return -1;
	}
	libWriter.write( m_libString.c_str(), m_libString.length() );
	libWriter.close();

	// calculate the coverage
	for( int i = 0; i < (int) myContigs->size(); i++ ){
		myContigs->at( i )->CalculateCov();
	}

	return 1;
}

// calculate mean and std using Interquartile Range
int ContigConverter::IQR( vector<double> *dis, double &mean, double &std, string libName, double &percentageOfOutliers, double &numberOfUsedPairs,
		       PetLibrary *currentLibrary )
{
	ofstream disWriter( (Configure::OUTPUT_FOLDER + "distanceOfMatePairs_" + libName).c_str(), ios::out );
	if( disWriter == NULL ){
		cout<<"ERROR: Cannot open "<<(Configure::OUTPUT_FOLDER + "distanceOfMatePairs_" + libName)<<" file"<<endl;
		return -1;
	}

	//cerr<<"distance size: "<<dis->size()<<endl;
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

	for( int i = 0; i < (int) dis->size(); i++ ){
		disStr += itos( dis->at( i ) ) + "\n";
		if( dis->at( i ) >= lowerbound && dis->at( i ) <= upperbound ){
			sum += dis->at( i );
			currentLibrary->AddOccuOfDistance( dis->at( i ) );
			//distributionResults += itos( dis->at( i ) ) + "\n";
			num++;
		}
	}
	
	currentLibrary->CalculatePossibility();

	disWriter.write( disStr.c_str(), disStr.length() );
	percentageOfOutliers = (1.0 - num/(double)dis->size()) * 100;
	numberOfUsedPairs = num;

	mean = sum / num;
	std = 0;
	// calculate std
	for( int i = 0; i < (int) dis->size(); i++ ){
		if( dis->at( i ) >= lowerbound && dis->at( i ) <= upperbound ){
			std += pow( dis->at( i ) - mean, 2 );
		}
	}

	std = std / num;
	std = pow( std, 0.5 );

	return 1;
}

// check if two lines represent a pair of reads
bool ContigConverter::IsPair( string firstRead, string secondRead ){
	string firstName = firstRead.substr( 0, firstRead.length() -2 );
	string secondName = secondRead.substr( 0, secondRead.length() - 2 );

	if( firstName == secondName ){
		return true;
	}
	else
		return false;
}


// analyze velvet file
// return -1 if failed
int ContigConverter::AnalyzeVelvet( ifstream *contigReader, vector<Contig*> *contigs ){
	string line;
	
	string name;
	double length = 0;
	double cov = 0;

	while( getline( *contigReader, line ) ){
		if( line.at( 0 ) == '>' ){
			// save previous contig
			if( length > 0 ){
				Contig *newContig = new Contig( name, length, cov );
				m_contigNameHashMap->insert( pair<const char*, int>( name.c_str(), contigs->size() ) );
				contigs->push_back( newContig );
			}

			// start new contig
			name = line.substr( 1, line.length() - 1 );
			int pos = line.find_last_of( "_" );
			cov = atof( line.substr( pos + 1, line.length() - pos - 1 ).c_str() );
			length = 0;;
		}
		else{
			length += line.length();
		}
	}

	// save last contig
	Contig *newContig = new Contig( name, length, cov );
	m_contigNameHashMap->insert( pair<const char*, int>( name.c_str(), contigs->size() ) );
	contigs->push_back( newContig );

	return 1;
}

// calculate the distance of two reads in the same contig
// including the read length
int ContigConverter::CalculateReadDisOnSameContig( vector<string> *firstColumn, vector<string> *secondColumn ){
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

// analyze soapdenovo file
int ContigConverter::AnalyzeSoapDenovo( ifstream *contigReader, vector<Contig*> *contigs ){
	string line;
	
	string name;
	double length = 0;
	double cov = 1;

	while( getline( *contigReader, line ) ){
		if( line.length() == 0 )
			continue;
			
		if( line.at( 0 ) == '>' ){
			// save previous contig
			if( length > 0 ){
				Contig *newContig = new Contig( name, length, cov );
				m_contigNameHashMap->insert( pair<const char*, int>( name.c_str(), contigs->size() ) );
				contigs->push_back( newContig );
			}

			// start new contig
			vector<string> *content = new vector<string>;
			Split( line, " ", content );
			name = line.substr( 1, content->at( 0 ).length() - 1 );
			cov = atof( content->at( 1 ).c_str() );
			length = 0;
			content->clear();
			delete content;
		}
		else{
			length += line.length();
		}
	}

	// save last contig
	Contig *newContig = new Contig( name, length, cov );
	m_contigNameHashMap->insert( pair<const char*, int>( name.c_str(), contigs->size() ) );
	contigs->push_back( newContig );

	return 1;
}


// analyze soapdenovo file
int ContigConverter::AnalyzeSoapDenovoContig( ifstream *contigReader, vector<Contig*> *contigs ){
	string line;
	
	string name;
	double length = 0;
	double cov = 1;

	while( getline( *contigReader, line ) ){
		if( line.length() == 0 )
			continue;
			
		if( line.at( 0 ) == '>' ){
			// save previous contig
			if( length > 0 ){
				Contig *newContig = new Contig( name, length, cov );
				m_contigNameHashMap->insert( pair<const char*, int>( name.c_str(), contigs->size() ) );
				contigs->push_back( newContig );
			}

			// start new contig
			vector<string> *content = new vector<string>;
			Split( line, " ", content );
			name = line.substr( 1, content->at( 0 ).length() - 1 );
			Split( line, "_", content );
			cov = atof( content->at( 1 ).c_str() );
			length = 0;
			content->clear();
			delete content;
		}
		else{
			length += line.length();
		}
	}

	// save last contig
	Contig *newContig = new Contig( name, length, cov );
	m_contigNameHashMap->insert( pair<const char*, int>( name.c_str(), contigs->size() ) );
	contigs->push_back( newContig );

	return 1;
}


// analyze opera scaffold file 
int ContigConverter::AnalyzeOperaScaf( ifstream *contigReader, vector<Contig*> *contigs ){
        string line;

        string name;
        double length = 0;
        double cov = 1;

        while( getline( *contigReader, line ) ){
                if( line.length() == 0 )
                        continue;

                if( line.at( 0 ) == '>' ){
                        // save previous contig                                                                                                                                   
                        if( length > 0 ){
                                Contig *newContig = new Contig( name, length, cov );
                                m_contigNameHashMap->insert( pair<const char*, int>( name.c_str(), contigs->size() ) );
                                contigs->push_back( newContig );
				//cerr<<name<<"\t"<<length<<"\t"<<cov<<endl;
                        }

                        // start new contig                      
                        vector<string> *content = new vector<string>;
                        Split( line, "\t", content );
                        name = line.substr( 1, content->at( 0 ).length() - 1 );
                        vector<string> tempContent;
			Split( content->at( 2 ), " ", &tempContent );
                        cov = atof( tempContent.at( 1 ).c_str() );
                        length = 0;
                        content->clear();
                        delete content;
                }
		else{
			length += line.length();
		}
	}

	// save last contig                                                                                                                                                      
	Contig *newContig = new Contig( name, length, cov );
	m_contigNameHashMap->insert( pair<const char*, int>( name.c_str(), contigs->size() ) );
	contigs->push_back( newContig );

	return 1;
}


// analyze normal fasta file
int  ContigConverter::AnalyzeNormalFasta( ifstream *contigReader, vector<Contig*> *contigs )
{
	string line;
	
	string name;
	double length = 0;
	double cov = 1;

	while( getline( *contigReader, line ) ){
		if( line == "" )
			continue;

		if( line.at( 0 ) == '>' ){
			// save previous contig
			if( length > 0 ){
				Contig *newContig = new Contig( name, length, cov );
				m_contigNameHashMap->insert( pair<const char*, int>( name.c_str(), contigs->size() ) );
				contigs->push_back( newContig );
				//cerr<<name<<endl;
			}

			// start new contig
			char nameChar[ 500 ];
			sscanf( line.c_str(), "%s", nameChar );
			name = string( nameChar );
			name = name.substr( 1, name.length() );
			length = 0;
		}
		else{
			length += line.length();
		}
	}

	// save last contig
	Contig *newContig = new Contig( name, length, cov );
	m_contigNameHashMap->insert( pair<const char*, int>( name.c_str(), contigs->size() ) );
	contigs->push_back( newContig );
	return 1;
}

// filter repeat
void ContigConverter::FilterRepeat( vector<Contig*> *contigs, Graph *graph ){
	vector<Contig*> usedContigs;		// save the non-repeat long contigs

	if( Configure::HAPLOID_COVERAGE == -1 ){
		// user did not specify the hyploid coverage, need to calculate
		double sum = 0;
		double num = 0;
		for( int i = 0; i < (int) contigs->size(); i++ ){
			if( contigs->at( i )->GetCov() > 0 && contigs->at( i )->GetLength() > Configure::CONTIG_SIZE_THERSHOLD ){
				sum += contigs->at( i )->GetCov() * contigs->at( i )->GetLength();
				num += contigs->at( i )->GetLength();
			}
		}

		double average = sum / num;

		Configure::HAPLOID_COVERAGE = average / Configure::PLOIDY;
	}

	double threshold = Configure::HAPLOID_COVERAGE * (Configure::PLOIDY + 0.5);

	//cerr<<"average coverage is "<<average<<endl;

#ifdef DEBUG
	cout<<"average coverage is "<<Configure::HAPLOID_COVERAGE<<endl;
#endif

	double totalUniqueContig = 0;
	double totalRepeat = 0;
	double totalSmallContig = 0;
	vector<Contig*>::iterator contigIter = contigs->begin();
	while( contigIter != contigs->end() ){
		(*contigIter)->SetIfRepeat( false );

		//Set the contig as repeat
		if( ! Configure::KEEP_REPEAT_FULL && (*contigIter)->GetCov() > threshold ){
			// remove repeat contigs
			m_repeatContigs->push_back( *contigIter );
			//(*contigIter)->SetContigType( REPEAT );
			(*contigIter)->SetIfRepeat( true );
			totalRepeat++;

			graph->InsertRepeatContigName( *contigIter );
		}
		
		/*
		if( (*contigIter)->GetLength() < Configure::SUPER_SMALL_CONTIG_THRESHOLD ) {
			// remove super small contigs
			m_smallContigs->push_back( *contigIter );
			(*contigIter)->SetContigType( SUPER_SMALL );
			totalSmallContig++;
		}
		else */
		if( (*contigIter)->GetLength() < Configure::CONTIG_SIZE_THERSHOLD ){
			// remove small contigs
			m_smallContigs->push_back( *contigIter );
			(*contigIter)->SetContigType( SMALL );
			totalSmallContig++;
		}
		else{
			//usedContigs.push_back( *contigIter );
			(*contigIter)->SetContigType( BIG );
			totalUniqueContig++;
		}

		usedContigs.push_back( *contigIter );
		contigIter++;
	}

	//cerr<<"Total unique contigs: "<<totalUniqueContig<<endl;
	//cerr<<"Total repeats: "<<totalRepeat<<endl;
	//cerr<<"Total small contigs: "<<totalSmallContig<<endl;

	// save to graph
	graph->SetContigNum( usedContigs.size() );
	for( int i = 0; i < (int) usedContigs.size(); i++ )
		graph->AddContig( usedContigs.at( i ) );
}

// analyze opera's contig file
int ContigConverter::AnalyzeStatistic( ifstream *contigReader, vector<Contig*> *contigs ){
	string line;

	vector<string> temp;
	getline( *contigReader, line );
	while( getline( *contigReader, line ) ){
		/*int pos = line.find( "\t" );
		if( pos == -1 ){
			return -1;
		}

		string name = line.substr( 0, pos );
		string length = line.substr( pos + 1, line.length() - pos - 1 );
		Contig *newContig = new Contig( name, atoi( length.c_str() ), 1 );
		contigs->push_back( newContig );*/

		Split( line, "\t", &temp );
		string name = temp.at( 0 );
		string length = temp.at( 1 );
		string cov = temp.at( 2 );
		Contig *newContig = new Contig( name, atoi( length.c_str() ), atoi(cov.c_str()) );
		contigs->push_back( newContig );

		temp.clear();
	}

	return 1;
}

// remove small contigs only
void ContigConverter::RemoveSmallContigs( vector<Contig*> *contigs, Graph *graph ){
	vector<Contig*> usedContigs;		// save the non-repeat long contigs

	vector<Contig*>::iterator contigIter = contigs->begin();
	while( contigIter != contigs->end() ){
		if( (*contigIter)->GetLength() < Configure::CONTIG_SIZE_THERSHOLD ){
			// remove small contigs
			m_smallContigs->push_back( *contigIter );
			(*contigIter)->SetContigType( SMALL );
		}
		else{
			//usedContigs.push_back( *contigIter );
			(*contigIter)->SetContigType( BIG );
		}

		usedContigs.push_back( *contigIter );
		contigIter++;
	}

	// save to graph
	graph->SetContigNum( usedContigs.size() );
	for( int i = 0; i < (int) usedContigs.size(); i++ )
		graph->AddContig( usedContigs.at( i ) );
}

// print list of contigs
int ContigConverter::PrintContigs( list<Contig*> *contigs, string fileName ){
	ofstream contigWriter( fileName.c_str() );

	if( contigWriter == NULL ){
		cout<<"ERROR: Cannot open "<<fileName<<" file"<<fileName<<endl;
		return -1;
	}
	
	string head = "Contig ID\tLength\tCoverage\n";
	contigWriter.write( head.c_str(), head.length() );

	list<Contig*>::const_iterator iter;
	for( iter = contigs->begin(); iter != contigs->end(); iter++ ){
		string result = (*iter)->GetName(); 
		char buffer[500];  
		//sprintf( buffer, "\t%.0f\t%.0f\n", (*iter)->GetLength(), (*iter)->GetCov() + 0.5 );
		sprintf( buffer, "\t%.0f\t%.0f\n", (*iter)->GetLength(), (*iter)->GetCov() );
		string tempString( buffer );
		result += tempString;
		contigWriter.write( result.c_str(), result.length() );
	}

	contigWriter.close();

	return 1;
}

