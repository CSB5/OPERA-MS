#include "CommonFunction.h"


// convert int to string
string itos( int i ){
	string result;
	stringstream buf;
	buf<<i;
	buf>>result;
	buf.clear();
	return result;
}

void Split( string line, string label, vector<string> *column ){
	column->clear();
	char * cstr, *p;
	cstr = new char [line.length()+1];
	strcpy (cstr, line.c_str());

	p=strtok (cstr, label.c_str() );
	while (p!=NULL)
	{
		column->push_back( string( p ) );
		p=strtok(NULL, label.c_str() );
	}

	delete[] cstr;  
}

void Split( string line, string label, list<string> *column ){
	column->clear();
	int pos = line.find( label );
	while( pos != -1 ){
		column->push_back( line.substr( 0, pos ) );
		line = line.substr( pos + 1, line.length() - pos - 1 );
		pos = line.find( label );
	}
	if( line != "" )
		column->push_back( line );
}


// get the opposite orientation 
int GetOppositeOri( int ori ){
	if( ori == PLUS )
		return MINUS;
	else
		return PLUS;
}

// print thousand comma format of a number
string PrintNumber( double num ){
	string result;
	char buffer [50];
	if( num > 1000000000 ){
		// gb
		//printf( "%.2f Gb", num/1000000000 );
  	sprintf (buffer, "%.2f Gb", num/1000000000 );
	}
	else if( num > 1000000 ){
		// mb
		//printf( "%.2f Mb", num/1000000 );
		sprintf (buffer, "%.2f Mb", num/1000000 );
	}
	else if( num > 1000 ){
		// kb
		//printf( "%.2f Kb", num/1000 );
		sprintf (buffer, "%.2f Kb", num/1000 );
	}
	else{
		// bp
		//printf( "%d bp", num );
		sprintf (buffer, "%d bp", (int)num );
	}
	
	result = string( buffer );
	return result;
}

// print thousand comma format of a number to file
void PrintNumberToFile( FILE *file, double num ){
	if( num > 1000000000 ){
		// gb
		fprintf( file, "%.2f Gb", num/1000000000 );
	}
	else if( num > 1000000 ){
		// mb
		fprintf( file, "%.2f Mb", num/1000000 );
	}
	else if( num > 1000 ){
		// kb
		fprintf( file, "%.2f Kb", num/1000 );
	}
	else{
		// bp
		fprintf( file, "%d bp", (int)num );
	}
}

// check if it is a number (coverage)
bool IsNumber( string s ){
	for( int i = 0; i < (int)s.size(); i++ ){
		char x = s.at( i );
		if( ( x >= '0' && x <= '9' ) || x == '.' )
			continue;
		else
			return false;
	}

	return true;
}

// convert sam format to bowtie format
void ConvertSamToBowtie( vector<string> *line )
{
	// Sam format: read id, ori, scaffold id, position, , , , , ,read sequence, mapping quality,
	// bowtie format: read id, ori, scaffold id, position, read sequence,

	line->at( 3 ) = itos( atoi( line->at( 3 ).c_str() ) - 1 );

	line->at( 4 ) = line->at( 9 );
	int sign = atoi( line->at( 1 ).c_str() );
	if( sign & 0X10 )
	{ // reverse sequence
		line->at( 1 ) = "-";
	}
	else
	{
		line->at( 1 ) = "+";
	}
}

// check if the file is sam format
bool IsSAM( string fileName )
{
	bool result = true;
	string tempFileName = fileName;
	if( fileName.substr( fileName.length() - 4, 4 ) == ".bam" ){
		//cerr<<"is sam\n";
		// it is a bam format
		tempFileName = Configure::OUTPUT_FOLDER + "temporary.sam";
		string cmd = "mkfifo " + tempFileName;
		if(system( cmd.c_str() ) );
		//cerr<<"finish mkfifo\n";
		cmd  = Configure::SAMDIR.c_str() + std::string("samtools view ") + fileName + " > " + tempFileName + "&"; 
		if(system( cmd.c_str() ) );
		//cerr<<"finish samtools\n";
	}

	//cerr<<"open file\n";
	//ifstream mapReader( fileName.c_str() );
	ifstream mapReader( tempFileName.c_str() );

	if( mapReader == NULL )
	{
		cout<<"error reading mapping file"<<endl;
		return false;
	}
	//cerr<<"succeed\n";
		
	string line;							
	getline( mapReader, line );
	
	vector<string> *temp = new vector<string>;
	Split( line, "\t", temp );
	//cerr<<"split line\n";
	if( temp->size() <= 1 || temp->at( 1 ) == "+" || temp->at( 1 ) == "-" )
		result = false;

	temp->clear();
	delete temp;

	if( fileName.substr( fileName.length() - 4, 4 ) == ".bam" ){
		string cmd = "rm " + Configure::OUTPUT_FOLDER + "temporary.sam";
		if(system( cmd.c_str() ) );
	}
	
	return result;
}


