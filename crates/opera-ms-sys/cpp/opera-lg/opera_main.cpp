#include "opera.h"

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
		exit(0);
		//return -1;
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
		//cerr<<"Samtools path "<<Configure::SAMDIR<<endl;
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
	//cerr << " *** " << m_opera->m_graph->GetContigIndex( string("k99_2"));
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
	//cerr << " *** " << m_opera->m_graph->GetContigIndex( string("k99_16"));
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
	// explicit: this main() is also called as a function by the Rust worker
	return 0;
}

