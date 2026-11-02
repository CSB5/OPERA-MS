#include "ARTreeNode.h"

ARTreeNode::ARTreeNode(void)
{
	mainMap = new map<int, ARTreeNode*>;
	edgeSet = new map<TreeNode*, int>;
}

ARTreeNode::ARTreeNode( int value ){
	edgeSet = new map<TreeNode*, int>;
	mainMap = new map<int, ARTreeNode*>;
	m_value = value;
}

ARTreeNode::~ARTreeNode(void)
{
	delete mainMap;
	delete edgeSet;
}

// check if dangling edge has existed
// if not, insert and return false
// else, return true
bool ARTreeNode::insertDE( TreeNode *node, int unhappyEdges ){
	map<TreeNode*, int>::iterator iter;
	if( edgeSet->empty() )
		iter = edgeSet->end();
	else
		iter = edgeSet->find( node );

	if( iter == edgeSet->end() ){
		// creat the connection
		edgeSet->insert( pair<TreeNode*, int>( node, unhappyEdges ) );
		return false;
	}
	else{
		if( (*iter).second <= unhappyEdges )
			return true;
		else{
			edgeSet->erase( iter );
			edgeSet->insert( edgeSet->end(), pair<TreeNode*, int>( node, unhappyEdges ) );
			return false;
		}
	}
}

// get a list of all treenode
list<ARTreeNode*>* ARTreeNode::GetChildren(){
	list<ARTreeNode*> *result = new list<ARTreeNode*>;
	for( map<int, ARTreeNode*>::iterator iter = mainMap->begin(); iter != mainMap->end(); iter++ ){
		result->push_back( (*iter).second );
	}
	return result;
}

// insert a value
// return the treenode which contain the key
ARTreeNode* ARTreeNode::insert( int key ){
	map<int, ARTreeNode*>::iterator iter = mainMap->find( key );
	if( iter == mainMap->end() ){
		// not exist
		ARTreeNode *newNode = new ARTreeNode( key );
		iter = mainMap->insert( iter, pair<int, ARTreeNode*>( key, newNode ) );
		return (*iter).second;
	}
	else
		return (*iter).second;
}
