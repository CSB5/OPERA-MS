#include "Tree.h"

Tree::Tree(void)
{
	arHead = new ARTreeNode( -1 );
	deHead = new TreeNode( -1 );
	m_totalNum = 0;
}

Tree::~Tree(void)
{
	delete arHead;
	delete deHead;
}


// check if the partial scaffold exists
// return true if exist, false if not
bool Tree::IfExist( string *arString, string *deString, int numOfUnhappyEdges ){
	// get active region value
	vector<string> *content = new vector<string>;
	Split( *arString, "\t", content );
	list<int> *arValue = new list<int>;
	for( int i = 0; i < (int)content->size() / 2 ; i++ ){
		int value = atoi( content->at( i * 2 ).c_str() ) + 1;		// avoid 0, add 1 to all values
		if( content->at( i * 2 + 1 ) == "-" )
			value = -value;
		arValue->push_back( value );
	}
	content->clear();
	// get last node of active region
	ARTreeNode *arNode = Insert( arValue, arHead );


	// get dangling edges value
	Split( *deString, "\t", content );
	list<int> *deValue = new list<int>;
	for( int i = 0; i < (int)content->size(); i++ ){
		deValue->push_back( atoi( content->at( i ).c_str() ) );
	}
	content->clear();
	// get last node of dangling edges
	TreeNode *deNode = Insert( deValue, deHead );

	delete content;
	delete arValue;
	delete deValue;

	// check if they are connected and the number of unhappy edges is smaller
	if( arNode->insertDE( deNode, numOfUnhappyEdges ) == false ){
		m_totalNum++;
		return false;
	}
	else{
		//cout<<"exist"<<endl;
		return true;
	}
}

// insert a new sequence, return corresponding last treenode
TreeNode* Tree::Insert( list<int> *sequence, TreeNode *head ){
	TreeNode *resultNode = head;
	for( list<int>::iterator iter = sequence->begin(); iter != sequence->end(); iter++ ){
		resultNode = resultNode->insert( *iter );
	}
	return resultNode;
}

// insert a new sequence, return corresponding last treenode
ARTreeNode* Tree::Insert( list<int> *sequence, ARTreeNode *head ){
	ARTreeNode *resultNode = head;
	for( list<int>::iterator iter = sequence->begin(); iter != sequence->end(); iter++ ){
		resultNode = resultNode->insert( *iter );
	}
	return resultNode;
}

// clear the tree
void Tree::Clear(){
	m_totalNum = 0;
	// clear active region
	list<ARTreeNode*> *tempARList = new list<ARTreeNode*>;
	tempARList->push_back( arHead );
	while( !tempARList->empty() ){
		ARTreeNode *arNode = tempARList->front();
		tempARList->pop_front();

		// save the children
		list<ARTreeNode*> *children = arNode->GetChildren();
		tempARList->insert( tempARList->end(), children->begin(), children->end() );
		delete children;
		
		delete arNode;
	}
	delete tempARList;

	// clear dangling edge
	list<TreeNode*> *tempDEList = new list<TreeNode*>;
	tempDEList->push_back( deHead );
	while( !tempDEList->empty() ){
		TreeNode *deNode = tempDEList->front();
		tempDEList->pop_front();

		// save the children
		list<TreeNode*> *children = deNode->GetChildren();
		tempDEList->insert( tempDEList->end(), children->begin(), children->end() );
		delete children;
		
		delete deNode;
	}
	delete tempDEList;
	
	// creat new heads
	arHead = new ARTreeNode( -1 );
	deHead = new TreeNode( -1 );
}

// get the total number of partial scaffolds
int Tree::GetTotalNum()
{
	return m_totalNum;
}
