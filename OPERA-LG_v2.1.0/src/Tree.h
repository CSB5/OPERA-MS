#pragma once

#include "TreeNode.h"
#include "ARTreeNode.h"
#include "CommonFunction.h"

#include <string>
#include <list>
#include <iostream>

using namespace std;

class Tree
{
public:
	Tree(void);
	~Tree(void);

	// attributes
public:
	ARTreeNode *arHead;			// the head of active region
	TreeNode *deHead;			// the head of active region, the first node is empty
	int m_totalNum;

	// methods
public:
	// check if the partial scaffold exists
	bool IfExist( string *arString, string *deString, int numOfUnhappyEdges );
	// clear the tree
	void Clear();
	// get the total number of partial scaffolds
	int GetTotalNum();

private:
	// insert a new sequence, return corresponding last treenode
	inline TreeNode* Insert( list<int> *sequence, TreeNode *head );
	// insert a new sequence, return corresponding last treenode
	inline ARTreeNode* Insert( list<int> *sequence, ARTreeNode *head );

};
