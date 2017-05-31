#pragma once
#include "TreeNode.h"
#include <set>

using namespace std;
class ARTreeNode
{
public:
	ARTreeNode(void);
	ARTreeNode( int value );
	~ARTreeNode(void);

	// attributes
public:
	map<int, ARTreeNode*> *mainMap;			// the mapping of the same class, such as active region or edges
	int m_value;
	map<TreeNode*, int> *edgeSet;			// the set of all possible dangling edges
												// record the number of unhappy edges
	// method
public:
	bool insertDE( TreeNode *node, int unhappyEdges );
	// get a list of all treenode
	list<ARTreeNode*>* GetChildren();
	// insert a value
	ARTreeNode* insert( int key );
};
