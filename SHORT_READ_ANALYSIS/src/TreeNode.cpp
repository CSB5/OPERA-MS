#include "TreeNode.h"

TreeNode::TreeNode(void)
{
	mainMap = new map<int, TreeNode*>;
}

TreeNode::TreeNode( int value ){
	mainMap = new map<int, TreeNode*>;
	m_value = value;
}

TreeNode::~TreeNode(void)
{
	delete mainMap;
}

// insert a value
// return the treenode which contain the key
TreeNode* TreeNode::insert( int key ){
	map<int, TreeNode*>::iterator iter = mainMap->find( key );
	if( iter == mainMap->end() ){
		// not exist
		TreeNode *newNode = new TreeNode( key );
		iter = mainMap->insert( iter, pair<int, TreeNode*>( key, newNode ) );
		return (*iter).second;
	}
	else
		return (*iter).second;
}

// get a list of all treenode
list<TreeNode*>* TreeNode::GetChildren(){
	list<TreeNode*> *result = new list<TreeNode*>;
	for( map<int, TreeNode*>::iterator iter = mainMap->begin(); iter != mainMap->end(); iter++ ){
		result->push_back( (*iter).second );
	}
	return result;
}
