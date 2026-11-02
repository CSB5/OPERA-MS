#pragma once

#include <map>
#include <list>


using namespace std;
class TreeNode
{
public:
	TreeNode(void);
	TreeNode( int value );
	~TreeNode(void);

	// properties
public:
	map<int, TreeNode*> *mainMap;			// the mapping of the same class, such as active region or edges
	int m_value;

	// methods
public:
	// insert a value
	TreeNode* insert( int key );
	// get a list of all treenode
	list<TreeNode*>* GetChildren();
};
